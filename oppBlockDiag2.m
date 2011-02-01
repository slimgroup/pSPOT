classdef oppBlockDiag2 < oppSpot
    %OPPBLOCKDIAG   Operator-diagonal operator in parallel sans overlap
    %
    %   B = oppBlockDiag(OP1, OP2,...,OPN,GATHER) creates a compound
    %   block operator with the input operators OP1, OP2,... on the diagonal of
    %   B, e.g., B = DIAG([OP1 OP2 ... OPN]). When multiplying the operators
    %   are distributed among the labs and multiplied locally on each.
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %
    %   B = opBlockDiag(WEIGHT,OP1,...,OPN,GATHER) additionally
    %   weights each block by the elements of the vector WEIGHT. If
    %   only a single operator is given it is replicated as many times
    %   as there are weights.
    %
    %   B = opBlockDiag(N,OP,GATHER) similar as above with WEIGHT
    %   equal to ones(N,1). This will cause operator OP to be repeated
    %   N times.
    %
    %   B = opBlockDiag(WEIGHT,A,GATHER) where A is a 3D numerical
    %   matrix. This will split the slices of A along the diagonal of B.
    %
    %   See also opFoG, opKron, opDictionary.
    
    %   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
    %   See the file COPYING.txt for full copyright information.
    %   Use the command 'spot.gpl' to locate this file.
    
    %   http://www.cs.ubc.ca/labs/scl/spot
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        weights;
        dist3D; % Distributed Explicit 3D matrix (assume distributed along 3rd dim)
        locm; % Local sizes for use when multiplying
        locn;
        chicodist; % Codistributor for the child operators
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppBlockDiag2(varargin)
            
            opdist3D = 0; % Temporary flag for working with distributed 3D
            opcflag = 0;
            localm = 0;
            localn = 0;
            % Extract gather
            if isscalar( varargin{end} ) && nargin > 1
                gather = varargin{end};
                varargin(end) = [];
            else
                gather = 0;
            end
            
            % Extract weight
            nargs = length(varargin);
            if ~isnumeric(varargin{1})
                weights = ones(nargs,1);
            else
                weights = varargin{1};
                if isempty(weights), weights = 1; end;
                [m,n] = size(weights);
                
                % Repeating Op
                if spot.utils.isposintscalar(weights)&& nargs == 2
                    weights = ones(weights,1);
                    varargin(1) = [];
                    
                    % Normal
                elseif (((m == 1) && (n == nargs-1)) || ...
                        ((n == 1) && (m == nargs-1)) || ...
                        ((m == 1) && (n == 1)))
                    weights = ones(nargs-1,1).*weights(:);
                    varargin(1) = [];
                    
                    % 3D input matrix
                elseif nargs == 1 && ndims(varargin{1}) == 3
                    nSlices = size(varargin{1},3);
                    if isdistributed(varargin{1}), opdist3D = 1; end
                    weights = ones(size(varargin{1},3),1);
                elseif nargs == 2 && ndims(varargin{2}) == 3
                    nSlices = size(varargin{2},3);
                    if isdistributed(varargin{2}), opdist3D = 1; end
                    if ((m == 1) && (n == nSlices)) || ...
                            (m == nSlices) && (n == 1)
                        weights = ones(nSlices,1).*weights(:);
                    else
                        weights = ones(nSlices,1);
                    end
                    varargin(1) = [];
                    
                else
                    weights = ones(nargs,1);
                end
            end
            
            % Extracting and distributing the operators (in distributed cell
            % arrays)
            if ~exist('nSlices','var')
                % Check for empty operators and remove them
                ops = cellfun(@isempty,varargin);
                assert(any(~ops),'At least one operator must be specified.');
                arrayfun(@(ind) warning('input "%d" is empty',ind), find(ops));
                varargin(ops) = [];
                if isdistributed(varargin) && iscell(varargin)
                    opList = varargin;
                else
                    opList = distributed(varargin);
                end
            else        %3D matrix case
                if opdist3D % Distributed 3D Matrix (kept as an array)
                    opList = varargin{1};
                else % Undistributed 3D Matrix (stuffed into distributed cell)
                    opList = distributed.cell( 1, nSlices);
                    for i = 1:nSlices
                        opList(i) = {varargin{1}(:,:,i)};
                    end
                end
            end
            clear varargin;
            
            % Check for pSpot operators
            if ~opdist3D % Distributed 3D has no spot operators
                for ops = cellfun(@(p) isa(p,'oppSpot'), opList)
                    if ops, error('oppSpot operators are not supported');end
                end
            end
            
            % Convert all Arguments to operators -except 3D matrix case
            if ~exist('nSlices','var')
                ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
                opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));
            end
            
            % Check complexity and setup single Op cases
            if exist('nSlices','var')           % 3D matrix
                if opdist3D % Distributed 3D matrix
                    [m,n,nSlices] = size(opList);
                    m = m*nSlices;  n = n*nSlices;
                    for i = 1:nSlices
                        if ~isreal(opList(:,:,i)), opcflag = 1; end
                    end
                else
                    [m,n] = cellfun( @size, opList ); %have to use cellfun because
                    % distributed cells don't support {} subscripting
                    m = m*nSlices;    n = n*nSlices;
                    for i = 1:length(opList)
                        if ~cellfun(@isreal,opList(i)), opcflag = 1; end
                    end
                end
                opcflag = opcflag || ...
                    ~all( isreal( weights ) );
                linear = 1;
                
            elseif length(opList) == 1    % Repeat 1 operator
                
                [m,n] = cellfun( @size, opList );
                m = m * length(weights);
                n = n * length(weights);
                opcflag  = ~cellfun( @isreal, opList) || ~all( isreal( weights));
                linear = cellfun(@(opl) opl.linear, opList );
                [opList{1:length(weights)}] = deal(opA);
                clear opA
                
            else
                
                [m,n] = cellfun(@size,opList);
                m = sum(m);     n = sum(n);
                real = cellfun(@isreal,opList);
                opcflag = ~all(real);
                linear = cellfun(@(p) logical(p.linear), opList);
                linear = all(linear);
                
            end
            
            % Convert distributed attributes to non-distributed scalars
            if isdistributed(m) || isdistributed(opcflag) ||...
                    isdistributed(linear) || opdist3D
                spmd
                    if iscodistributed(m), m = getLocalPart(m); end
                    if iscodistributed(n), n = getLocalPart(n); end
                    if iscodistributed(opcflag), opcflag = getLocalPart(opcflag); end
                    if iscodistributed(linear), linear = getLocalPart(linear); end
                    
                    % Extract local m and n for future use
                    if opdist3D % Distributed 3D case
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
                end
                clear childs; clear mm; clear nn; clear ss; clear i;
                
                % Getting the scalar out of the composites
                if isa(m,'Composite'), m = m{1}; end; m = m(1);
                if isa(n,'Composite'), n = n{1}; end; n = n(1);
                if isa(opcflag,'Composite'), opcflag = opcflag{1}; end
                if isa(linear,'Composite'), linear = linear{1}; end
                
            end % Conversion
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', m, n);
            op.locm = localm;
            op.locn = localn;
            op.cflag    = opcflag;
            op.linear   = linear;
            op.children = opList;
            op.weights  = weights;
            op.sweepflag= true;
            op.gather   = gather;
            op.dist3D = opdist3D;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = 'pBlockDiag(';
            childs = op.children;
            if op.dist3D % For explicit distributed 3D matrices
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
                    if ~cellfun(@isnumeric, childs(1) )
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
            
            % Multidimensional x preprocessing
            if size(x,2) > 1, error('Please vectorize your x!'); end
            numdims = size(x,1)/op.n;
            if numdims > 1
                for i = 0:numdims-1 % Loop through every dimension on x
                    y = op.multiply(x((op.n*i)+1:op.n*(i+1)),mode); % Recursive calls
                end
                return;
            end
            
            % Get the number of labs running
            nlabs = findResource;
            nlabs = nlabs.ClusterSize;
            
            % Setting up partition variables
            partition = [];
            finpartition = []; % partition for the final answer
            
            % Setting up class variables
            opweights = op.weights;
            childs = op.children;
            opdist3D = op.dist3D;
            opgather = op.gather;
            
            % Forward mode (mode = 1)
            if mode == 1
                % Distributing x into the correct chunks
                % Process multidimensional x (assuming x is always passed in
                % vectorized)
                glosize = [op.n 1]; % Setting up the global sizes
                finsize = [op.m 1];
                for i = 1:nlabs % summing local n for partition
                    partition = [partition sum(op.locn{i})];
                    % Setting up the partition for the answers
                    finpartition = [finpartition sum(op.locm{i})];
                end
                
            else % mode == 2
                % Distributing x into the correct chunks
                % Process multidimensional x (assuming x is always passed in
                % vectorized)
                glosize = [op.m 1];
                finsize = [op.n 1];
                for i = 1:nlabs
                    partition = [partition sum(op.locm{i})];
                    finpartition = [finpartition sum(op.locn{i})];
                end
            end % end of mode 2
            
            % Check for distribution of x and redistribute if necessary
            % Also distributes weights
            if isdistributed(x)
                spmd
                    defcodist = getCodistributor(x);
                    codist = codistributor1d(1,partition,glosize);
                    if defcodist.Partition ~= codist.Partition;
                        warning('WARNING: x is not distributed correctly!');
                        x = redistribute(x,codist);
                    end
                end
            else    % Distribute x
                warning('WARNING: x is not distributed!');
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
                childs = getLocalPart(childs);
                x = getLocalPart(x);
                if opdist3D % multiply for distributed 3D matrix case
                    for i = 1:size(childs,3) % Extracting operators and spotify
                        opchilds{i} = opMatrix(childs(:,:,i));
                    end
                    ind = globalIndices(codist,3);
                else
                    ind = globalIndices(codist);
                end
                B = opBlockDiag(opweights(ind),opchilds{:});
                if mode == 1
                    tmpy = B*x;
                else
                    tmpy = B'*x;
                end
                
                if opgather % Gather all data to lab 1
                    fincodist = codistributor1d(1,finpartition,finsize);
                    tmpy = codistributed.build(tmpy,fincodist);
                    tmpy = gather(tmpy,1); % Slow and crude version -- for now
                end
                
            end % spmd
            if op.gather
                y = tmpy{1}; %gathered data is on lab 1
            else
                y = tmpy;
            end % gather
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef