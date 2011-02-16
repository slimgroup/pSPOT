classdef oppNumBlockDiag < oppSpot
    %OPPBLOCKDIAG   Operator-diagonal operator in parallel sans overlap
    %   Supports distributed x and distributed 3D Matrices
    %
    %   B = oppNumBlockDiag(OP1, OP2,...,OPN,GATHER) creates a compound
    %   block operator with the input operators OP1, OP2,... on the diagonal of
    %   B, e.g., B = DIAG([OP1 OP2 ... OPN]). When multiplying the operators
    %   are distributed among the labs and multiplied locally on each.
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %
    %   B = oppNumBlockDiag(WEIGHT,OP1,...,OPN,GATHER) additionally
    %   weights each block by the elements of the vector WEIGHT. If
    %   only a single operator is given it is replicated as many times
    %   as there are weights.
    %
    %   B = oppNumBlockDiag(N,OP,GATHER) similar as above with WEIGHT
    %   equal to ones(N,1), where N is a positive integer. This will cause
    %   operator OP to be repeated N times.
    %
    %   B = oppNumBlockDiag([WEIGHT],A,GATHER) where A is a 3D numerical
    %   matrix. This will slice A along the 3rd dimension and use the 2D
    %   slices as the blocks for the block-diagonal operator B.
    %   Note: 2D x vector not supported for numeric 3D
    %
    %   B = oppNumBlockDiag([WEIGHT],A,NUMCOLS_X,GATHER) where A is a 3D numerical
    %   matrix. This will slice A along the 3rd dimension and use the 2D
    %   slices Ak to build Kronecker operator kron(opDirac(NUMCOLS_X), opMatrix(Ak)).
    %   These will then be used as the blocks for the block-diagonal operator B.
    %   In this mode the flag GATHER is REQUIRED to distinguish this from the
    %   previous case.
    %
    %
    %   See also oppBlockDiag, oppDictionary, oppStack
    
    %   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
    %   See the file COPYING.txt for full copyright information.
    %   Use the command 'spot.gpl' to locate this file.
    
    %   http://www.cs.ubc.ca/labs/scl/spot
    
    % Tim's change Feb 1, 2011: Now works correctly when number of workers exceed number of operators!! Damn you thomas T_T I spent 2 hours
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        locm; % Local sizes for use when multiplying
        locn;
        numcols; % number of implied RHS when handling x that is implicitly multidimensional, default to 1
        isBlockMatrix;  % a flag that indicates whether the block diagonal operators are built from slices of a 3D array
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppNumBlockDiag(varargin)
            
            % Check Matlabpool
            if matlabpool('size') == 0
                error('Matlabpool is not on');
            end
            
            % Setting up the variables
            localm = 0;
            localn = 0;
            Nrep = [];
            numcols_x = 1;
            gather = 0;
            isBlockMatrix = false;
            
            % Arguments Extraction
            % Extract gather & numcols
            if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot') % Gather
                if isscalar(varargin{end-1}) && ~isa(varargin{end},'opSpot') % numcols
                    numcols_x = varargin{end-1};
                    varargin(end-1) = [];
                end
                gather = varargin{end};
                varargin(end) = [];
            end
            
            % Extract weight
            nargs = length(varargin);
            if ~isnumeric(varargin{1}) || ndims(varargin{1})==3   % No weights
                weights = ones(nargs,1);
            else
                if isscalar(varargin{1}) 
                    if ndims(varargin{2}) ==3
                        weights = varargin{1};
                    else% N for repeating ops
                        Nrep = varargin{1};
                    end
                else                       % weights
                    weights = varargin{1};
                    if isempty(weights), weights = 1; end;
                end
                varargin(1) = []; % Remove
                nargs = nargs - 1;
            end
            
            if nargs == 1 % 3D input matrix or repeating ops
                if ndims(varargin{1}) == 3 % 3D Matrix
                    nSlices = size(varargin{1},3);
                    Cube = varargin{1};
                    isBlockMatrix = true;
                    if isscalar(weights)
                        weights = weights.*ones(size(varargin{1},3),1); % default
                    else
                        if length(weights(:)) == nSlices
                            weights = ones(nSlices,1).*weights(:);
                        else
                            error('WEIGHTS size mismatch');
                        end
                    end
                elseif ~isempty(Nrep) % Repeating Op
                    weights = ones(Nrep,1);
                else
                    % Single spot operator repeated weights times
                    weights = ones(length(weights(:)),1).*weights(:);
                end
                
            else % Normal
                if (((length(weights(:)) == nargs)) || length(weights(:)) == 1)
                    weights = ones(nargs,1).*weights(:);
                else
                    error('WEIGHTS size mismatch');
                end
            end
            
            % Operators Extraction and Processing
            % Extracting and distributing the operators (in distributed cell
            % arrays)
            if ~isBlockMatrix % Normal or repeating ops
                % Check for empty operators and remove them
                ops = cellfun(@isempty,varargin);
                assert(any(~ops),'At least one operator must be specified.');
                arrayfun(@(ind) warning('No:Input','input "%d" is empty',ind), find(ops));
                varargin(ops) = [];
                opList = distributed(varargin); % There is no way varargin is distributed
                
            else % 3D matrix case
                if isdistributed(Cube) % Distributed 3D Matrix
                    warn3D = {0};
                    spmd
                        numcodist = getCodistributor(Cube);
                        if numcodist.Dimension ~= 3
                            warn3D = 1;
                            Cube = redistribute(Cube, codistributor1d(3));
                        end
                    end
                    if warn3D{1}
                        warning('WarnDist:Wrongdistribution',...
                                    '3D Matrix is not distributed correctly!\nCommencing automatic redistribution...');
                    end
                    clear warn3D;
                    opList = Cube;
                else % Undistributed 3D Matrix
                    warning('WarnDist:Nodistribution',...
                        '3D Matrix is not distributed!\nCommencing automatic redistribution...');
                    opList = distributed(Cube);
                end
            end
            clear varargin;
            
            % Check for pSpot operators
            if ~isBlockMatrix % Numeric matrix has no spot operators
                for ops = cellfun(@(p) isa(p,'oppSpot'), opList)
                    if ops, error('oppSpot operators are not supported');end
                end
            end
            
            % Convert all Arguments to operators -except 3D matrix case
            if ~isBlockMatrix
                ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
                if isdistributed(ops)
                    spmd
                        ops = getLocalPart(ops);
                    end
                    ops = ops{:};
                end
                opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));
            end
            
            % Check complexity and setup single Op cases
            if isBlockMatrix  % 3D matrix
                
                [m,n,nSlices] = size(opList);
                m = m*nSlices;  n = n*nSlices;
                cflag = ~isreal(opList) || ~all(isreal(weights));
                linear = 1;
                % numcols_x for multiple RHS
                m = m*numcols_x;    n = n*numcols_x;
                
            elseif length(opList) == 1    % Repeat 1 operator
                
                [m,n] = cellfun( @size, opList );
                m = m * length(weights);
                n = n * length(weights);
                cflag  = ~cellfun( @isreal, opList) || ~all( isreal( weights));
                linear = cellfun(@(opl) opl.linear, opList );
                spmd, opList = getLocalPart(opList); end; % De-distributify opList
                opList = opList{1};
                for i = 1:length(weights)
                    opList{i} = opList{1};
                end
                opList = distributed(opList);
            else
                
                [m,n] = cellfun(@size,opList);
                m = sum(m);     n = sum(n);
                real = cellfun(@isreal,opList);
                cflag = ~all(real);
                linear = cellfun(@(p) logical(p.linear), opList);
                linear = all(linear);
                
            end
            
            % Post-processing and Construction of Operator
            % Convert distributed attributes to non-distributed scalars
            spmd
                % Extract local m and n for future use
                if isBlockMatrix % Distributed 3D case
                    childs = getLocalPart(opList);
                    [mm,nn,ss] = size(childs);
                    for i = 1:ss
                        localm(1,i) = mm;
                        localn(1,i) = nn;
                    end
                else
                    childs = getLocalPart(opList);
                    [localm,localn] = cellfun(@size,childs);
                end
                if iscodistributed(m), m = getLocalPart(m); end
                if iscodistributed(n), n = getLocalPart(n); end
                if iscodistributed(cflag), cflag = getLocalPart(cflag); end
                if iscodistributed(linear), linear = getLocalPart(linear); end
            end
            clear childs; clear mm; clear nn; clear ss; clear i;
            % Getting the scalar out of the composites
            if isa(m,'Composite'), m = m{1}; end; m = m(1);
            if isa(n,'Composite'), n = n{1}; end; n = n(1);
            if isa(cflag,'Composite'), cflag = cflag{1}; end
            if isa(linear,'Composite'), linear = linear{1}; end
            
            % Construct operator
            op = op@oppSpot('pnumBlockDiag', m, n);
            op.locm = localm;
            op.locn = localn;
            op.cflag    = cflag;
            op.sweepflag = false;
            op.linear   = linear;
            op.children = opList;
            op.weights  = weights;
            op.sweepflag= true;
            op.gather   = gather;
            op.isBlockMatrix = isBlockMatrix;
            if exist('numcols_x','var')
                op.numcols = numcols_x;
            else
                op.numcols = 1;
            end
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = 'numBlockDiag(';
            childs = op.children;
            if op.isBlockMatrix % For explicit distributed 3D matrices
                spmd
                    childs = getLocalPart(childs);
                    string = '';
                    [m,n,nSlices] = size(childs);
                    for i=1:nSlices
                        string = [string,'Matrix(',int2str(m),',',int2str(n),'), '];
                    end
                end % spmd
            else
                spmd
                    childs = getLocalPart(childs);
                    string = '';
                    if ~isempty(childs)
                        if ~cellfun(@isnumeric, childs(1))
                            for childs = childs
                                string = [string,char(childs{1}),', '];
                            end
                        else
                            [m,n] = cellfun(@size, childs);
                            for i=1:length(childs)
                                string = [string,'Matrix(',int2str(m(i)),',', ...
                                    int2str(n(i)),'), '];
                            end
                        end
                    end
                end % spmd
            end % if dist3D
            
            % Concatenating composite string
            stringy = '';
            for i = 1:length(string)
                stringy = [stringy string{i}];
            end
            str = [str,stringy(1:end-2),')'];
        end % Display
        
    end % Methods
    
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
                                 
            % Setting up partition variables
            partition = [];
            finpartition = []; % partition for the final answer
            
            % Setting up class variables
            opweights = op.weights;
            childs = op.children;
            opgather = op.gather;
            isBlockMatrix = op.isBlockMatrix;
            numcols = op.numcols;
            
            % Matrix x preprocessing
            if isBlockMatrix
                if size(x,2) ~= 1
                    error('Please vectorize your x');
                else
                    multicols = 1;
                end
            else
                multicols = size(x,2); % for matrix multiplications
            end
            
            
            % Get the number of labs running
            if matlabpool('size') == 0
                warning('Now:LookieHere','Matlabpool is not open');
                nlabs = 1;
            else
                nlabs = matlabpool('size');
            end
            
            % Forward mode (mode = 1)
            if mode == 1
                % Distributing x into the correct chunks
                glosize = [op.n multicols]; % Setting up the global sizes
                finsize = [op.m multicols];
                for i = 1:nlabs % summing local n for partition
                    if i <= size(op.locn,2)
                        partition = [partition sum(op.locn{i})];
                        % Setting up the partition for the answers
                        finpartition = [finpartition sum(op.locm{i})];
                    else
                        partition = [partition 0];
                        finpartition = [finpartition 0];
                    end
                end
                
            else % mode == 2
                glosize = [op.m multicols];
                finsize = [op.n multicols];
                for i = 1:nlabs
                    if i <= size(op.locm,2)
                        partition = [partition sum(op.locm{i})];
                        finpartition = [finpartition sum(op.locn{i})];
                    else
                        partition = [partition 0];
                        finpartition = [finpartition 0];
                    end
                end
            end % end of mode 2
            
            % take into account implied multiple RHS in multidimensional x
            partition = partition .* op.numcols;
            finpartition = finpartition .* op.numcols;
            
            % Check for distribution of x and redistribute if necessary
            if isdistributed(x)
                warnx = {0};
                spmd
                    defcodist = getCodistributor(x);
                    codist = codistributor1d(1,partition,glosize);
                    if defcodist.Partition ~= codist.Partition
                        warnx = 1;
                        x = redistribute(x,codist);
                    end
                end
                if warnx{1}
                    warning('WarnDist:Wrongdistribution',...
                                'x is not distributed correctly!\nCommencing automatic redistribution...');
                end
                clear warnx;
            else    % Distribute x
                warning('WarnDist:Nodistribution',...
                    'x is not distributed!\nCommencing automatic redistribution...');
                spmd
                    codist = codistributor1d(1,partition,glosize);
                    x = codistributed(x,codist);
                end
            end
            clear codist; % codistributor for x cleared
            clear partition; % partition for x cleared
            
            % Multiplication starts
            spmd
                
                % Multiply using opBlockDiag locally
                codist = getCodistributor(childs);
                opchilds = getLocalPart(childs);
                local_x = getLocalPart(x);
                
                if isBlockMatrix % multiply for distributed 3D matrix case
                    ind = globalIndices(codist,3);
                else
                    ind = globalIndices(codist,2);
                end
                
                if ~isempty(opchilds)
                    
                    if isBlockMatrix
                        
                        num_child_ops = size(opchilds,3); % number of blocks local to this machine
                        local_weights = opweights(ind);
                        % reshape x if we have multiple RHS
                        % assertion: this block only executes for 3D matrix input. Therefore, all the sizes (m and n) of all
                        % child operators must be the same at this point.
                        local_x = reshape(local_x, [], numcols, num_child_ops);
                        if mode == 1
                            tmpy = zeros(size(opchilds,1), numcols, num_child_ops);
                            for k = 1:num_child_ops;
                                A = opchilds(:,:,k);
                                tmpy(:,:,k) = local_weights(k) .* (A * local_x(:,:,k));
                            end
                        else
                            tmpy = zeros(size(opchilds,2), numcols, num_child_ops);
                            for k = 1:num_child_ops;
                                A = opchilds(:,:,k);
                                tmpy(:,:,k) = conj(local_weights(k)) .* (A' * local_x(:,:,k));
                            end
                        end
                        
                        tmpy = tmpy(:);
                        
                    else
                        B = opBlockDiag(opweights(ind),opchilds{:});
                        if mode == 1
                            tmpy = B*local_x;
                        else
                            tmpy = B'*local_x;
                        end
                        
                    end
                    
                else
                    % local part of distributed vector for nodes that do not contain data is of size (0,1) for some reason
                    tmpy = zeros(0,multicols);
                end
                
                fincodist = codistributor1d(1,finpartition,finsize);
                tmpy = codistributed.build(tmpy,fincodist,'noCommunication');
                
            end % spmd
            if mode == 1
                if op.gather == 1 || op.gather == 2
                    y = gather(tmpy);
                else
                    y = tmpy;
                end
            else % mode == 2
                if op.gather == 1 || op.gather == 3
                    y = gather(tmpy);
                else
                    y = tmpy;
                end
            end % gather
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef