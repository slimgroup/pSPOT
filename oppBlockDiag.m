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
            
            % Standard checking and setup sizes
            [opList,m,n,cflag,linear] = stdpspotchk(varargin{:});
            m = sum(m);     n = sum(n);
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', m, n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
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
            if ~isnumeric( op.children{1} )
                for child = op.children
                    str = [str,char(child{1}),', '];
                end
            else
                [m,n] = cellfun(@size, op.children);
                for i=1:length(op.children)
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
            opchildren = distributed(op.children);
            opweights = op.weights;
            gather = op.gather;
            opm = op.m;
            opn = op.n;
            
            % Checking for distributed x
            if ~isdistributed(x)
                error('X is not distributed');
            end
            
            spmd
                % Setting up the local parts
                codist = getCodistributor(opchildren);
                wind = globalIndices(codist,2); % local weights indices
                local_weights = opweights(wind);
                local_children = getLocalPart(opchildren);
                local_x = getLocalPart(x);
                
                % Setting up codistributor for y
                finpart = codistributed.zeros(1,numlabs);
                if mode ==  1
                    fingsize = [opm size(x,2)]; % final partition
                else
                    fingsize = [opn size(x,2)]; % final global size
                end
                
                if ~isempty(local_children)
                    if length(local_children) == 1 % Just to avoid repeating op case
                        B = local_children{1};
                        if mode == 1
                            y = local_weights .* (B*local_x);
                            finpart(labindex) = B.m;
                        else
                            y = conj(local_weights) .* (B'*local_x);
                            finpart(labindex) = B.n;
                        end
                    else
                        B = opBlockDiag(local_weights,local_children{:});
                        if mode == 1
                            y = B*local_x;
                            finpart(labindex) = B.m;
                        else
                            y = B'*local_x;
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = divide(op,x,mode)
            % Setting up the variables
            opchildren = distributed(op.children);
            opweights = op.weights;
            gather = op.gather;
            opm = op.m;
            opn = op.n;
            
            % Checking for distributed x
            if ~isdistributed(x)
                error('X is not distributed');
            end
            
            spmd
                % Setting up the local parts
                codist = getCodistributor(opchildren);
                wind = globalIndices(codist,2); % local weights indices
                local_weights = opweights(wind);
                local_children = getLocalPart(opchildren);
                local_x = getLocalPart(x);
                
                % Setting up codistributor for y
                finpart = codistributed.zeros(1,numlabs);
                if mode ==  1
                    fingsize = [opn size(x,2)];
                else
                    fingsize = [opm size(x,2)];
                    local_weights = conj(local_weights);
                end
                
                if ~isempty(local_children)
                    
                    % Extracting local sizes
                    localm = 0; localn = 0;
                    for i=1:length(local_children)
                        child = local_children{i};
                        localm = localm + child.m;
                        localn = localn + child.n;
                    end
                    
                    if mode == 1
                        finpart(labindex) = localn;
                        
                        % Preallocate y
                        y = zeros(localn,size(x,2));
                        
                        % Divide
                        j = 0;
                        k = 0;
                        for i=1:length(local_children) % Divide by operator
                            child = local_children{i};
                            y(j+1:j+child.n,:) = local_weights(i) .* ...
                                (child \ local_x(k+1:k+child.m,:));
                            j = j + child.n;
                            k = k + child.m;
                        end
                        
                    else % mode 2
                        finpart(labindex) = localm;
                        
                        % Preallocate y
                        y = zeros(localm,size(x,2));
                        
                        % Divide
                        j = 0;
                        k = 0;
                        for i=1:length(local_children) % Divide by operator
                            child = local_children{i};
                            child = child';
                            y(j+1:j+child.n,:) = local_weights(i) .* ...
                                (child \ local_x(k+1:k+child.m,:));
                            j = j + child.n;
                            k = k + child.m;
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
            
        end % divide
    end % Protected Methods
    
end % Classdef









































