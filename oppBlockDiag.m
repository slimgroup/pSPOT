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
    %   B = oppBlockDiag([WEIGHT],OP1,...,OPN,GATHER) additionally
    %   weights each block by the elements of the vector WEIGHT. If
    %   only a single operator is given it is replicated as many times
    %   as there are weights. If a scalar weight is given for multiple
    %   operators the scalar WEIGHT will be expanded to the size of the
    %   operators, ie. WEIGHT*ones(N,1)
    %
    %   B = oppBlockDiag(N,OP,GATHER) similar as above with WEIGHT
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
            
            % Check Matlabpool
            assert(matlabpool('size') > 0, 'Matlabpool is not on');
            
            % Extract gather
            if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [];
            else
                gather = 0;
            end
            
            % Extract weights and fill in repeating ops
            nargs = length(varargin);
            
            if isnumeric(varargin{1}) % weights
                weights = varargin{1};
                weights = weights(:);
                
                if nargs == 2 % Repeating ops
                    % repeating N times
                    if spot.utils.isposintscalar(varargin{1}) 
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
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', sum(m), sum(n));
            op.cflag       = cflag;
            op.linear      = linear;
            op.children    = opList;
            op.weights     = weights;
            op.sweepflag   = true;
            op.gather      = gather;
            op.opsn        = n;
            op.opsm        = m;
            
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
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Drandn
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function x = drandn(A,Ncols)
%             ncols = 1;
%             if nargin == 2 % for easy multivectoring
%                 ncols = Ncols;
%             end
%             
%             % Distribute children
%             opchildren = distributed(A.children);
%             
%             spmd, chicodist = getCodistributor(opchildren); end
%             
%             chicodist = chicodist{1};
%             chipart = chicodist.Partition;
%             chinum = 0;
%             for i=1:matlabpool('size')
%                 xpart(i) = 0;
%                 for j=chinum+1:chinum+chipart(i)
%                     child = A.children{j};
%                     xpart(i) = xpart(i) + child.n;
%                 end
%                 chinum = chinum + chipart(i);
%             end
%             xgsize = [A.n ncols];
%             
%             n = A.n;
%             if isreal(A)
%                 spmd
%                     xcodist = codistributor1d(1,xpart,xgsize);
%                     x = codistributed.randn(n,ncols,codistributor1d(1));
%                     x = redistribute(x,xcodist);
%                 end
%             else
%                 spmd
%                     xcodist = codistributor1d(1,xpart,xgsize);
%                     x = codistributed.randn(n,ncols,codistributor1d(1)) +...
%                         1i*codistributed.randn(n,ncols,codistributor1d(1));
%                     x = redistribute(x,xcodist);
%                 end
%             end
%             
%         end % drandn
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Rrandn
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function x = rrandn(A,Ncols)
%             ncols = 1;
%             if nargin == 2 % for easy multivectoring
%                 ncols = Ncols;
%             end
%             
%             opchildren = distributed(A.children);
%                         
%             spmd, chicodist = getCodistributor(opchildren); end
%             
%             chicodist = chicodist{1};
%             chipart = chicodist.Partition;
%             chinum = 0;
%             for i=1:matlabpool('size')
%                 xpart(i) = 0;
%                 for j=chinum+1:chinum+chipart(i)
%                     child = A.children{j};
%                     xpart(i) = xpart(i) + child.m;
%                 end
%                 chinum = chinum + chipart(i);
%             end
%             xgsize = [A.m ncols];
%             
%             m = A.m;
%             
%             if isreal(A)
%                 spmd
%                     xcodist = codistributor1d(1,xpart,xgsize);
%                     x = codistributed.randn(m,ncols,codistributor1d(1));
%                     x = redistribute(x,xcodist);
%                 end
%             else
%                 spmd
%                     xcodist = codistributor1d(1,xpart,xgsize);
%                     x = codistributed.randn(m,ncols,codistributor1d(1)) +...
%                         1i*codistributed.randn(m,ncols,codistributor1d(1));
%                     x = redistribute(x,xcodist);
%                 end
%             end
%             
%         end % rrandn
        
    end % Methods
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Checking x
            assert(isdistributed(x),'x is not distributed');
            
            % Checking size of x
            spmd, xcodist   = getCodistributor(x); end
            xcodist   = xcodist{1};
            xpart     = xcodist.Partition;
            chipart   = pSPOT.utils.defaultDistribution(length(op.children));
            nlabs     = matlabpool('size');
            
            assert(xcodist.Dimension == 1,... % Dimensional check
                'x is not distributed along dimension 1');
            
            chinum = 0;
            for i=1:nlabs
                childm = sum(op.opsm(chinum+1:(chinum+chipart(i))));
                childn = sum(op.opsn(chinum+1:(chinum+chipart(i))));
                
                if mode == 1
                    assert(childn == xpart(i),...
                        'x size mismatch at lab %d, check distribution',i);
                else % mode 2
                    assert(childm == xpart(i),...
                        'x size mismatch at lab %d, check distribution',i);
                end
                chinum = chinum + chipart(i);
            end            
            
            % Setting up the variables and partition size
            loc_children = pSPOT.utils.compositeDef(op.children);
            loc_weights  = pSPOT.utils.compositeDef(op.weights);            
            if mode ==  1
                fingsize = [op.m size(x,2)]; % final global size
            else
                fingsize = [op.n size(x,2)];
            end
            
            spmd
                % Setting up the local parts
                loc_x    = getLocalPart(x);
                
                % Setting up codistributor for y
                finpart  = codistributed.zeros(1,numlabs);
                
                if ~isempty(loc_children)                    
                    if length(loc_children) == 1
                        % Just to avoid repeating op case
                        B = loc_children{1};
                        if mode == 1
                            y = loc_weights .* (B*loc_x);
                            finpart(labindex) = B.m;
                        else
                            y = conj(loc_weights) .* (B'*loc_x);
                            finpart(labindex) = B.n;
                        end
                    else
                        B = opBlockDiag(loc_weights,loc_children{:});
                        if mode == 1
                            y = B*loc_x;
                            finpart(labindex) = B.m;
                        else
                            y = B'*loc_x;
                            finpart(labindex) = B.n;
                        end
                        
                    end
                else
                    y = zeros(0,size(x,2));
                end
                
                % Check for sparsity
                aresparse            = codistributed.zeros(1,numlabs);
                aresparse(labindex)  = issparse(y);
                if any(aresparse), y = sparse(y); end;
                
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
            gather = op.gather;
            
            % Checking x
            assert(isdistributed(x),'x is not distributed');
            
            % Checking size of x
            spmd, xcodist   = getCodistributor(x); end
            xcodist   = xcodist{1};
            xpart     = xcodist.Partition;
            chipart   = pSPOT.utils.defaultDistribution(length(op.children));
            nlabs     = matlabpool('size');
            
            assert(xcodist.Dimension == 1,... % Dimensional check
                'x is not distributed along dimension 1');
            
            chinum = 0;
            for i=1:nlabs
                childm = sum(op.opsm(chinum+1:(chinum+chipart(i))));
                childn = sum(op.opsn(chinum+1:(chinum+chipart(i))));
                
                if mode == 1
                    assert(childm == xpart(i),...
                        'x size mismatch at lab %d, check distribution',i);
                else % mode 2
                    assert(childn == xpart(i),...
                        'x size mismatch at lab %d, check distribution',i);
                end
                chinum = chinum + chipart(i);
            end
            
            % Setting up the variables and global sizes            
            loc_children = pSPOT.utils.compositeDef(op.children);            
            if mode ==  1
                fingsize = [op.n size(x,2)];
                loc_weights  = pSPOT.utils.compositeDef(op.weights);
            else
                fingsize = [op.m size(x,2)];
                loc_weights = pSPOT.utils.compositeDef(conj(op.weights));
            end
            
            spmd
                % Setting up the local parts
                loc_x = getLocalPart(x);
                
                % Setting up codistributor for y
                finpart = codistributed.zeros(1,numlabs);
                
                if ~isempty(loc_children)                    
                    % Extracting local sizes
                    sizemn = cellfun(@size,loc_children,'UniformOutput',false);
                    sizemn = [sizemn{:}];
                    localm = sum(sizemn(1:2:end)); % sum odd
                    localn = sum(sizemn(2:2:end)); % sum even
                    
                    if mode == 1
                        finpart(labindex) = localn;
                        
                        % Preallocate y
                        y = zeros(localn,size(x,2));
                        
                        % Divide
                        j = 0; k = 0;
                        for i=1:length(loc_children) % Divide by operator
                            child = loc_children{i};
                            y(j+1:j+child.n,:) = loc_weights(i) .* ...
                                (child \ loc_x(k+1:k+child.m,:));
                            j = j + child.n;
                            k = k + child.m;
                        end
                        
                    else % mode 2
                        finpart(labindex) = localm;
                        
                        % Preallocate y
                        y = zeros(localm,size(x,2));
                        
                        % Divide
                        j = 0; k = 0;
                        for i=1:length(loc_children) % Divide by operator
                            child = ctranspose(loc_children{i});
                            y(j+1:j+child.n,:) = loc_weights(i) .* ...
                                (child \ loc_x(k+1:k+child.m,:));
                            j = j + child.n;
                            k = k + child.m;
                        end                        
                    end
                else
                    y = zeros(0,size(x,2));
                end
                
                % Check for sparsity
                aresparse           = codistributed.zeros(1,numlabs);
                aresparse(labindex) = issparse(y);
                if any(aresparse),y = sparse(y); end;
                
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