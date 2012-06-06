classdef oppDictionary < oppSpot
    %OPPDICTIONARY  Dictionary of concatenated operators that acts in
    %               parallel
    %
    %   D = oppDictionary([WEIGHTS],OP1,OP2,...OPn,GATHER) creates a
    %   dictionary operator consisting of the concatenation of all
    %   operators. Each of the operators is multiplied with the
    %   corresponding part of x on seperate labs.
    %
    %   D = oppDictionary(N,OP) creates a stacked operator A using N number
    %   of repeating operators OP.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results
    %
    %   Optional WEIGHTS vector:
    %
    %               [WEIGHT1*OP1
    %                WEIGHT2*OP2
    %                   ...
    %                WEIGHTn*OPn]
    %
    %   If the same weight is to be applied to each operator, set
    %   WEIGHTS to a scalar. When WEIGHTS is empty [], it is set to
    %   one. The WEIGHT parameter can be omitted as long as OP1 is a
    %   Spot operator; if not there is no way to
    %   decide whether it is a weight vector or operator.
    %
    %   *Note - only spot operators can be used in an oppDictionary, pSpot
    %   operators act in parallel already, and cannot be used with
    %   oppDictionary.
    %
    %   **Note - As of now the operators will be distributed according to
    %   the Matlab default codistribution scheme. The x vector should be
    %   distributed with this in mind.
    %   for more info, type 'help codistributor1d'
    %
    %   See also oppBlockDiag, oppNumBlockDiag, oppStack
    
    %   Nameet Kumar - Oct 2010
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDictionary(varargin)
            
            % Check Matlabpool
            if matlabpool('size') == 0
                error('Matlabpool is not on');
            end
            
            % Setting up the variables
            gather = 0;
            
            % Check for gather parameter
            if isscalar( varargin{end} ) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [            ];
            end
            
            % Check for weights
            nargs = length(varargin);
            
            if isnumeric(varargin{1}) % weights
                weights = varargin{1};
                weights = weights(:);
                
                if nargs == 2 % Repeating ops
                    
                    if spot.utils.isposintscalar(varargin{1}) % repeating N times
                        weights = ones(weights,1);                        
                    end % Else: Repeating as many times as there are weights
                    
                    for i = 3:length(weights)+1
                        varargin{i} = varargin{2};
                    end
                    
                else % Non-repeating ops
                    
                    if isscalar(varargin{1}) % Same weight applied to all
                        weights = weights*ones(nargs-1,1);
                        
                    else
                        if length(varargin{1}) ~= nargs-1
                            % Incorrect weight size
                            error('Weights size mismatch');
                        end
                        % Else: Normal weights with normal ops
                    end
                end
                varargin(1) = []; % delete weights
                
            else    % no weights
                weights = ones(nargs,1);
            end
            
            % Standard pSpot checking and setup sizes
            [opList,m,n,cflag,linear] = ...
                pSPOT.utils.stdpspotchk(varargin{:});
            assert( all(m == m(1)), 'Operator sizes are not consistent');
            
            % Construct
            op = op@oppSpot('pDictionary', m(1), sum(n));
            op.cflag       = cflag;
            op.linear      = linear;
            op.children    = opList;
            op.weights     = weights;
            op.sweepflag   = true;
            op.gather      = gather;
            op.precedence  = 1;
            op.ddistscheme = n;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = ['[',char(op.children{1})];
            
            for ops=op.children(2:end)
                str = [str, ', ', char(ops{1})];
            end
            
            str = [str, ']'];
            
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = double(op)
            %OPPDICTIONARY.DOUBLE Distributed doubling of oppDictionary
            %    A = double(op) will apply double to each child operator
            %    in oppDictionary, and return a distributed dictionary
            %    of explicit operators.
            
            % Find the default partition
            chidist = pSPOT.utils            .defaultDistribution(length(op.children));
            chicomp = Composite();
            glosize = [op.m op.n];
            ind     = 1;
            
            % Send operators to corresponding labs
            for i = 1:length(chicomp)
                chicomp{i} = op.children(ind:(ind-1+chidist(i)));
                ind = ind + chidist(i);
            end
            
            % Actual doubling and rebuilding in spmd
            spmd
                childn          = 0;
                partition       = codistributed.zeros(1,numlabs);
                A               = zeros(glosize(1),childn);
                if ~isempty(chicomp)
                    for i = 1:length(chicomp)
                        child   = chicomp{i};
                        childn  = childn + child.n;
                    end
                    A           = zeros(glosize(1),childn);
                    partition(labindex) = childn;
                    k           = 0;
                    for i = 1:length(chicomp)
                        child   = chicomp{i};
                        n       = child.n;
                        A(:,k+1:k+n) = double(child);
                        k       = k+n;
                    end
                end
                partition       = gather(partition);
                codist = codistributor1d(2,partition,glosize);
                A = codistributed.build(A,codist,'noCommunication');
            end % spmd
            
        end % double
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = drandn(A,Ncols)
            if nargin == 2 % for easy multivectoring
                ncols = Ncols;
            else
                ncols = 1;
            end
            
            % Forward mode
            chipart   = pSPOT.utils.defaultDistribution(length(A.children));
            childnum  = 0;
            xpart     = zeros(1,matlabpool('size'));
            for i=1:matlabpool('size')
                xpart(i) = 0;
                for j=childnum+1:childnum+chipart(i)
                    child    = A.children{j};
                    xpart(i) = xpart(i) + child.n;
                end
                childnum     = childnum + chipart(i);
            end
            xgsize           = [A.n ncols];
            
            if isreal(A)
                spmd
                    xcodist = codistributor1d(1,xpart,xgsize);
                    x = codistributed.randn(xgsize,codistributor1d(1));
                    x = redistribute(x,xcodist);
                end
            else
                spmd
                    xcodist = codistributor1d(1,xpart,xgsize);
                    x = codistributed.randn(xgsize,codistributor1d(1)) +...
                        1i*codistributed.randn(xgsize,codistributor1d(1));
                    x = redistribute(x,xcodist);
                end
            end
            
        end % drandn
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rrandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = rrandn(A,Ncols)
            ncols = 1;
            if nargin == 2 % for easy multivectoring
                ncols = Ncols;
            end
            
            m = A.m;
            if isreal(A)
                x = randn(m,ncols);
            else
                x = randn(m,ncols) + 1i*randn(m,ncols);
            end            
        end % rrandn
        
    end % Methods    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 2 % Use oppStack, since a transpose of dictionary is
                % equivalent to a stack with transposed operators
                tchild = cellfun(@ctranspose,op.children,...
                    'UniformOutput', false);
                
                B = oppStack(opEye(op.n,op.m)); % Pseudo copy constructor
                B.children    = tchild;
                B.cflag       = op.cflag;
                B.sweepflag   = op.sweepflag;
                B.linear      = op.linear;
                B.gather      = op.gather;
                B.weights     = conj(op.weights); % Conj for complex numbers
                B.rdistscheme = cellfun(@(x) size(x,1),tchild);
                
                % Multiply
                if isdistributed(x)
                    tmpx = gather(x);
                else
                    tmpx = x;
                end
                y = B*tmpx;
                clear B;
                clear tmpx;
                return;
            end % Mode 2
            
            % X must be distributed
            assert(isdistributed(x),'X must be distributed');
            
            % Checking size of x
            spmd, xcodist = getCodistributor(x); end
            xcodist = xcodist{1};
            xpart   = xcodist.Partition;
            chipart = pSPOT.utils.defaultDistribution(length(op.children));
            
            assert(xcodist.Dimension == 1,... % Dimensional check
                'x is not distributed along dimension 1');
            
            % Size checkings
            childnum = 0;
            for i=1:matlabpool('size') 
                childn = sum(op.ddistscheme(childnum+1:(childnum+chipart(i))));
                assert(childn == xpart(i),...
                    'x size mismatch at lab %d, check your distribution',i);
                childnum       = childnum + chipart(i);
            end
            opchildren  = pSPOT.utils.compositeDef(op.children);
            loc_weights = pSPOT.utils.compositeDef(op.weights);
            
            % Mode 1
            % Setting up preallocation size
            ysize = [op.m size(x,2)];
            % This "renaming" is required to avoid passing in the whole op,
            % which for some weird reason stalls spmd
            
            spmd
                % Setting up local parts
                local_x = getLocalPart(x);
                                
                % Preallocate y
                y = zeros(ysize);
                
                if ~isempty(opchildren)
                    for i=1:length(opchildren)
                        opchildren{i} = loc_weights(i) * opchildren{i};
                    end
                    B = opDictionary(opchildren{:});
                    y = B*local_x;
                end
                
                % Check for sparsity
                aresparse            = codistributed.zeros(1,numlabs);
                aresparse(labindex)  = issparse(y);
                if any(aresparse), y = sparse(y); end;
                
                % Summing the results and distribute
                y = pSPOT.utils.global_sum(y); % The result now sits on lab 1
                y = codistributed(y,1,codistributor1d(2));
                
            end %spmd
            
            if op.gather
                y = gather(y);
            end % if we gathered, the data is on master client
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef