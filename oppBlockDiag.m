classdef oppBlockDiag < oppSpot
    %OPPBLOCKDIAG   Operator-diagonal operator in parallel.
    %
    %   B = oppBlockDiag(OP1, OP2,...,OPN,GATHER) creates a compound
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
    %   B = oppBlockDiag([WEIGHT],OP1,...,OPN,OVERLAP,GATHER) additionally
    %   weights each block by the elements of the vector WEIGHT. If
    %   only a single operator is given it is replicated as many times
    %   as there are weights. If a scalar weight is given for multiple
    %   operators the scalar WEIGHT will be expanded to the size of the
    %   operators, ie. WEIGHT*ones(N,1)
    %
    %   B = oppBlockDiag(N,OP,OVERLAP,GATHER) similar as above with WEIGHT
    %   equal to ones(N,1). This will cause operator OP to be repeated
    %   N times.
    %
    %   See also oppNumBlockDiag, oppDictionary, oppStack
    
    %   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
    %   See the file COPYING.txt for full copyright information.
    %   Use the command 'spot.gpl' to locate this file.
    
    %   http://www.cs.ubc.ca/labs/scl/spot
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        weights;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppBlockDiag(varargin)
            
            % Settin' up the variables
            gather = 0;
            
            % Extract gather
            if isscalar(varargin{end})
                gather = varargin{end};
                varargin(end) = [];
            end
            
            % Extract weights and fill in repeating ops
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
            
            % Check for empty operators and remove them
            ops = ~cellfun(@isempty,varargin);
            assert(any(ops),'At least one operator must be specified.');
            arrayfun(@(ind) warning('input "%d" is empty',ind), find(~ops));
            opList = varargin(ops);
            
            % Check for pSpot operators
            ops = cellfun(@(p) isa(p,'oppSpot'), opList);
            assert(~any(ops),' oppSpot operators are not supported');
            
            % Convert all non-Spot arguments to operators
            ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
            opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));
            
            % Check complexity and setup sizes
            [m,n] = cellfun(@size,opList);
            m = sum(m);     n = sum(n);
            real = cellfun(@isreal,opList);
            cflag = ~all(real);
            linear = cellfun(@(p) logical(p.linear), opList);
            linear = all(linear);
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', m, n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = distributed(opList);
            op.weights  = weights;
            op.sweepflag= true;
            op.gather   = gather;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = 'pBlockDiag(';
            opchildren = gather(op.children);
            if ~isnumeric( opchildren{1} )
                for child = opchildren
                    str = [str,char(child{1}),', '];
                end
            else
                [m,n] = cellfun(@size, opchildren);
                for i=1:length(opchildren)
                    str = [str,'Matrix(',int2str(m(i)),',', ...
                        int2str(n(i)),'), '];
                end
            end
            str = [str(1:end-2), ')'];
        end % Display
        
    end % Methods
    
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Setting up the variables
            opchildren = op.children;
            opweights = op.weights;
            gather = op.gather;
            opm = op.m;
            opn = op.n;
            
            % Checking for distributed x
            if ~isdistributed(x)
                error('X is not distributed');
            end
            
            % Preallocate y
            y = zeros(op.m,size(x,2));
                        
            spmd
                % Setting up the local parts
                codist = getCodistributor(opchildren);
                wind = globalIndices(codist,2); % local weights indices
                local_weights = opweights(wind);
                childs = getLocalPart(opchildren);
                x = getLocalPart(x);
                
                % Setting up codistributor for y
                finpart = codistributed.zeros(1,numlabs);
                if mode ==  1
                    fingsize = [opm size(x,2)];
                else
                    fingsize = [opn size(x,2)];
                end
                
                if ~isempty(childs)
                    if length(childs) == 1 % Just to avoid repeating op case
                        B = childs{1};
                        if mode == 1
                            y = local_weights .* (B*x);
                            finpart(labindex) = B.m;
                        else
                            y = local_weights .* (B'*x);
                            finpart(labindex) = B.n;
                        end
                    else
                        B = opBlockDiag(local_weights,childs{:});
                        if mode == 1
                            y = B*x;
                            finpart(labindex) = B.m;
                        else
                            y = B'*x;
                            finpart(labindex) = B.n;
                        end
                        
                    end
                else
                    y = zeros(0,size(x,2));
                end
                
                % Codistribute y
                fincodist = codistributor1d(1,finpart,fingsize);
                y = codistributed.build(y,fincodist,'noCommunication');
                
            end % spmd
            
            % Gather
            if mode == 1
                if op.gather == 1 || op.gather == 2
                    y = gather(y);
                end
            else % mode == 2
                if op.gather == 1 || op.gather == 3
                    y = gather(y);
                end
            end % gather
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef














