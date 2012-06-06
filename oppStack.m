classdef oppStack < oppSpot
    %OPPSTACK  Stack of vertically concatenated operators in parallel
    %
    %   S = oppStack([WEIGHTS], OP1, OP2, ...OPn,GATHER) creates a stacked 
    %   operator A consisting of the vertical concatenation of all 
    %   operators. When applied the operators are divided amongst the labs 
    %   and applied locally on each lab.
    %
    %   S = oppStack(N,OP) creates a stacked operator A using N number of
    %   repeating operators OP.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results.
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
    %   *Note - only spot operators can be used in an oppStack, pSpot
    %   operators act in parallel already, and cannot be used with
    %   oppDictionary.
    %
    %   **Note - As of now the operators will be distributed according to
    %   the Matlab default codistribution scheme.
    %   for more info, type 'help codistributor1d'
    %
    %   See also oppBlockDiag, oppNumBlockDiag, oppDictionary
    
    %   Nameet Kumar - Oct 2010
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppStack(varargin)
            
            % Check Matlabpool
            if matlabpool('size') == 0
                error('Matlabpool is not on');
            end
            
            % Setting up the variables
            gather = 0;
            
            % Check for gather parameter
            if isscalar( varargin{end} ) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [];
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
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            assert( all(n == n(1)), 'Operator sizes are not consistant');
            
            % Construct
            op = op@oppSpot('pStack', sum(m), n(1));
            op.cflag       = cflag;
            op.linear      = linear;
            op.children    = opList;
            op.weights     = weights;
            op.sweepflag   = true;
            op.gather      = gather;
            op.precedence  = 1;
            op.rdistscheme = m;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = ['[',char(op.children{1})];
            
            for ops=op.children(2:end)
                str = [str, '; ', char(ops{1})];
            end
            
            str = [str, ']'];
            
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = double(op)
            %OPPSTACK.DOUBLE Distributed doubling of oppStack
            %    A = double(op) will apply double to each child operator
            %    in oppStack, and return a distributed stack of explicit
            %    operators.
            
            opchildren = pSPOT.utils.compositeDef(op.children);
            globalsize = [op.m op.n];
            spmd
                % Setup distribution stuffs
                childm          = 0;
                partition       = codistributed.zeros(1,numlabs);
                
                % Preallocate
                A               = zeros(childm,globalsize(2));
                if ~isempty(opchildren)
                    A           = double(opStack(opchildren{:}));
                    partition(labindex) = childm;
                end
                partition = gather(partition);
                codist = codistributor1d(1,partition,globalsize);
                A = codistributed.build(A,codist,'noCommunication');
            end % spmd
            
        end % double
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = drandn(A,Ncols)
            ncols = 1;
            if nargin == 2 % for easy multivectoring
                ncols = Ncols;
            end
            
            n = A.n;
            if isreal(A)
                x = randn(n,ncols);
            else
                x = randn(n,ncols) + 1i*randn(n,ncols);
            end
            
        end % drandn
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rrandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = rrandn(A,Ncols)
            % for easy multivectoring
            if nargin == 2 
                ncols = Ncols;
            else
                ncols = 1;
            end
            
            chipart = pSPOT.utils.defGlobInd(length(A.children));
            xpart   = zeros(1,matlabpool('size'));
            for i=1:matlabpool('size')
                xpart(i) = sum(A.rdistscheme([chipart{i}]));
            end
            xgsize = [A.m ncols];
            
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
            
        end % rrandn
        
    end % Methods    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 2 % Use oppDictionary, since a transpose of stack is
                % transpose operators
                tchild = cellfun(@ctranspose,op.children,...
                    'UniformOutput', false);
                
                B = oppDictionary(opEye(op.n,op.m)); % Pseudo copy constructor
                B.children    = tchild;
                B.cflag       = op.cflag;
                B.sweepflag   = op.sweepflag;
                B.linear      = op.linear;
                B.gather      = op.gather;
                B.weights     = conj(op.weights); % Conj for complex numbers                
                B.ddistscheme = cellfun(@(x) size(x,2),tchild);
                
                % Multiply
                y = B*x;
                clear B;
                return;
            end % Mode 2
            
            % Checking distribution of x
            if isdistributed(x)
                spmd, xcodist = getCodistributor(x); end
                xcodist = xcodist{1};
                assert(xcodist.Dimension == 1,...
                    'x cannot be distributed along first dimension');
            end
            
            % Mode 1
            % Setting up class variables and partition sizes
            loc_children = pSPOT.utils.compositeDef(op.children);
            loc_weights  = pSPOT.utils.compositeDef(op.weights);
            fingsize     = [op.m size(x,2)]; % final global size
            chibind      = pSPOT.utils.defGlobInd(length(op.children));
            finpart      = distributed.zeros(1,matlabpool('size'));
            
            for i=1:matlabpool('size')
                finpart(i) = sum(op.rdistscheme([chibind{i}]));
            end
            
            spmd
                % Preallocate y
                y = zeros(getLocalPart(finpart),fingsize(2));
                
                if ~isempty(loc_children)                    
                    % Multiply
                    for i=1:length(loc_children)
                        loc_children{i} = loc_weights(i) * loc_children{i};
                    end
                    B = opStack(loc_children{:});
                    y = B*x;
                else
                    y = zeros(0,fingsize(2));
                end
                
                % Check for sparsity
                aresparse = codistributed.zeros(1,numlabs);
                aresparse(labindex) = issparse(y);                
                if any(aresparse), y = sparse(y); end;
                
                % Concatenating the results and distribute
                finpart   = gather(finpart);
                fincodist = codistributor1d(1,finpart,fingsize);
                y = codistributed.build(y,fincodist,'noCommunication');                
            end %spmd
            
            if op.gather
                y = gather(y);
            end    %if we gathered, the data is on master client
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef