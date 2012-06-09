classdef oppNumBlockDiag < oppSpot
%OPPNUMBLOCKDIAG   Operator-diagonal operator in parallel sans overlap
%   Supports distributed x and distributed 3D Matrices
%
%   B = oppNumBlockDiag([WEIGHT],A,GATHER) where A is a 3D numerical
%   matrix. This will slice A along the 3rd dimension and use the 2D
%   slices as the blocks for the block-diagonal operator B.
%   Note: 2D x vector not supported for numeric 3D
%   GATHER specifies whether to gather the results to a local array
%   or leave them distributed, default is 0.
%   GATHER = 0 will leave them distributed.
%   GATHER = 1 will gather the results of forwards or adjoint 
%              multiplication.
%   GATHER = 2 will gather only in forward mode.
%   GATHER = 3 will gather only in backward (adjoint) mode.
%
%   B = oppNumBlockDiag([WEIGHT],A,ncols_x,GATHER) where A is a 3D numerical
%   matrix. This will slice A along the 3rd dimension and use the 2D
%   slices Ak to build Kronecker operator kron(opDirac(ncols_x), opMatrix(Ak)).
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
    ncols; % number of implied RHS when handling x that is implicitly multidimensional, default to 1

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
        assert(matlabpool('size') > 0, 'Matlabpool is not open');

        % Setting up the variables
        ncols_x = 1;
        gather  = 0;

        % Arguments Extraction
        % Extract gather & ncols
        if isscalar(varargin{end}) % Gather
            if isscalar(varargin{end-1}) % ncols
                ncols_x         = varargin{end-1};
                varargin(end-1) = [];
            end
            gather        = varargin{end};
            varargin(end) = [];
        end

        % Extract weights
        nargs = length(varargin);
        if isdistributed(varargin{1}) || ~isnumeric(varargin{1}) ||...
           ndims(varargin{1})==3   % No weights
            weights = ones(nargs,1);
        else
            % weights
            weights = varargin{1};
            if isempty(weights), weights = 1; end;
            varargin(1) = []; % Remove
        end

        % Check weights
        nSlices = size(varargin{1},3);
        Cube    = double(varargin{1});
        if isscalar(weights)
            weights = weights.*ones(size(varargin{1},3),1); % default
        else
            if length(weights(:)) == nSlices
                weights = ones(nSlices,1).*weights(:);
            else
                error('WEIGHTS size mismatch');
            end
        end

        % Operators Extraction and Processing
        % Check for distribution of the 3D Matrix
        assert(isdistributed(Cube),'Cube must be distributed!');
        clear varargin;

        % Check complexity 
        [m,n,nSlices]     = size(Cube);
        m = m*nSlices;  n = n*nSlices;
        cflag  = ~isreal(Cube) || ~all(isreal(weights));
        linear = 1;
        % ncols_x for multiple RHS
        m = m*ncols_x;    n = n*ncols_x;


        % Post-processing and Construction of Operator
        % Convert distributed attributes to non-distributed scalars
        spmd
            % Check distribution dimension
            numcodist = getCodistributor(Cube);

            % Extract local m and n for future use
            childs     = getLocalPart(Cube);
            [mm,nn,ss] = size(childs);
            localm     = mm*ones(1,ss);
            localn     = nn*ones(1,ss);

            childs = []; mm = []; nn = []; ss = [];
        end
        numcodist = numcodist{1};
            assert(numcodist.Dimension == 3,...
                '3D Matrix is not distributed correctly!');

        % Construct operator
        op = op@oppSpot('pnumBlockDiag', m, n);
        op.locm      = localm;
        op.locn      = localn;
        op.cflag     = cflag;
        op.sweepflag = false;
        op.linear    = linear;
        op.children  = Cube;
        op.weights   = weights;
        op.sweepflag = true;
        op.gather    = gather;
        if exist('ncols_x','var')
            op.ncols = ncols_x;
        else
            op.ncols = 1;
        end
        op.opsn = (n/nSlices)*ones(1,nSlices);
        op.opsm = (m/nSlices)*ones(1,nSlices);

    end %Constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function str = char(op)
        % Initialize
        str = 'numBlockDiag(';
         childs = op.children;
%             if op.isBlockMatrix % For explicit distributed 3D matrices
            spmd
                childs = getLocalPart(childs);
                string = '';
                [m,n,nSlices] = size(childs);
                for i=1:nSlices
                    string = [string,'Matrix(',int2str(m),',',int2str(n),'), '];
                end
            end % spmd

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
        ncols = op.ncols;

        % Matrix x preprocessing
        if size(x,2) ~= 1
            error('Please vectorize your x');
        else
            multicols = 1;
        end


        % Get the number of labs running
        nlabs = matlabpool('size');

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
        partition = partition .* op.ncols;
        finpartition = finpartition .* op.ncols;

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
                error('x is not distributed correctly!');
            end
            clear warnx;
        else    % Distribute x
            error('x is not distributed!');
        end
        clear codist; % codistributor for x cleared
        clear partition; % partition for x cleared

        % Multiplication starts
        spmd

            % Multiply using opBlockDiag locally
            codist = getCodistributor(childs);
            opchilds = getLocalPart(childs);
            local_x = getLocalPart(x);
            ind = globalIndices(codist,3);

            if ~isempty(opchilds)

                num_child_ops = size(opchilds,3); % number of blocks local to this machine
                local_weights = opweights(ind);
                % reshape x if we have multiple RHS
                % assertion: this block only executes for 3D matrix input. Therefore, all the sizes (m and n) of all
                % child operators must be the same at this point.
                local_x = reshape(local_x, [], ncols, num_child_ops);
                if mode == 1
                    tmpy = zeros(size(opchilds,1), ncols, num_child_ops);
                    for k = 1:num_child_ops;
                        A = opchilds(:,:,k);
                        tmpy(:,:,k) = local_weights(k) .* (A * local_x(:,:,k));
                    end
                else
                    tmpy = zeros(size(opchilds,2), ncols, num_child_ops);
                    for k = 1:num_child_ops;
                        A = opchilds(:,:,k);
                        tmpy(:,:,k) = conj(local_weights(k)) .* (A' * local_x(:,:,k));
                    end
                end

                tmpy = tmpy(:);

            else
                % local part of distributed vector for nodes that do not contain data is of size (0,1) for some reason
                tmpy = zeros(0,multicols);
            end

            % Check for sparsity
            aresparse = codistributed.zeros(1,numlabs);
            aresparse(labindex) = issparse(tmpy);
            % labBarrier;
            if any(aresparse), tmpy = sparse(tmpy); end;

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