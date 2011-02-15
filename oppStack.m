classdef oppStack < oppSpot
    %OPPSTACK  Stack of vertically concatenated operators in parallel
    %
    %   oppStack([WEIGHTS], OP1, OP2, ...OPn,GATHER) creates a stacked operator
    %   consisting of the vertical concatenation of all operators. When applied
    %   the operators are divided amongst the labs and applied locally on each
    %   lab. 
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
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
    %   See also oppDictionary, opFoG, opSum.
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
                
                if isempty(weights) % Empty weights (why in the world???)
                    weights = 1;
                end
                
                if isscalar(varargin{1}) % Same weight applied to all
                    weights = weights*ones(nargs-1,1);
                    
                else
                    if length(varargin{1}) ~= nargs-1
                        % Incorrect weight size
                        error('Weights size mismatch');
                    end
                    % Else: Normal weights with normal ops
                end
                varargin(1) = []; % delete weights
                
            else    % no weights
                weights = ones(nargs,1);
            end
            
            % Standard pSpot checking and setup sizes
            [opList,m,n,cflag,linear] = stdpspotchk(varargin{:});
            assert( all(n == n(1)), 'Operator sizes are not consistant');
            m = sum(m);
            
            % Construct
            op = op@oppSpot('pStack', m, n(1));
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.weights  = weights;
            op.sweepflag= true;
            op.gather   = gather;
            op.precedence= 1;
            
            
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
            opchildren = distributed(op.children);
            childn = op.n;
            opm = op.m;
            spmd
                local_children = getLocalPart(opchildren);
                childm = 0;
                partition = codistributed.zeros(1,numlabs);
                globalsize = [opm childn];
                A = zeros(childm,childn);
                if ~isempty(local_children)
                    for i = 1:length(local_children)
                        child = local_children{i};
                        childm = childm + child.m;
                    end
                    A = zeros(childm,childn);
                    partition(labindex) = childm;
                    k = 0;
                    for i = 1:length(local_children)
                        child = local_children{i};
                        m = child.m;
                        A(k+1:k+m,:) = double(child);
                        k = k+m;
                    end
                end
                partition = gather(partition);
                codist = codistributor1d(1,partition,globalsize);
                A = codistributed.build(A,codist,'noCommunication');
            end % spmd
            
        end % double
        
    end % Methods
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 2 % Use oppDictionary, since a transpose of stack is
                
                opchildren = op.children; % equivalent to a dictionary with
                tchild = cell(1,length(opchildren)); % transposed operators
                for i = 1:length(opchildren) 
                    child = opchildren{i}; 
                    tchild{i} = child';
                end
                
                B = oppDictionary(opEye(op.n,op.m)); % Pseudo copy constructor
                B.children = tchild;
                B.cflag = op.cflag;
                B.sweepflag = op.sweepflag;
                B.linear = op.linear;
                B.gather = op.gather;
                B.weights = conj(op.weights); % Conj for complex numbers
                
                % Multiply
                y = B*x;
                clear B;
                return;
            end % Mode 2
            
            % Checking distribution of x
            if isdistributed(x)
                spmd, xcodist = getCodistributor(x); end
                xcodist = xcodist{1};
                if xcodist.Dimension == 1 % Checking distribution of x
                    error('x should not be distributed along first dimension');
                end
            end
            
            % Mode 1
            % Setting up class variables
            opchildren = distributed(op.children); % This "renaming" is
            opm = op.m; opn = op.n;   % required to avoid
            opweights = op.weights;   % passing in the whole op, which for
                                      % some weird reason stalls spmd
            spmd
                % Setting up local parts
                local_children = getLocalPart(opchildren);
                finpart = codistributed.zeros(1,numlabs); % final partition
                fingsize = [opm size(x,2)]; % final global size
                
                % Setting up weights
                codist = getCodistributor(opchildren);
                wind = globalIndices(codist,2);
                local_weights = opweights(wind);                
                                
                if ~isempty(local_children)
                    % Setup partition size
                    localm = 0;
                    for i=1:length(local_children)
                        child = local_children{i};
                        localm = localm + child.m;
                    end
                    finpart(labindex) = localm;
                    
                    % Preallocate y
                    y = zeros(localm,size(x,2));
                    
                    % Multiply
                    B = opStack(local_weights,local_children{:});
                    y = B*x;
                else
                    y = zeros(0,size(x,2));
                end
                
                % Concatenating the results and distribute
                finpart = gather(finpart);
                fincodist = codistributor1d(1,finpart,fingsize);
                y = codistributed.build(y,fincodist);
                
            end %spmd
            
            if op.gather
                y = gather(y); 
            end    %if we gathered, the data is on master client
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef