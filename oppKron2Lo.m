classdef oppKron2Lo < opKron
%oppKron2Lo  kronecker tensor product to act on a distributed vector. 
%
%   oppKron2Lo(A,B, gather)A and B are Spot operators or numeric matrices.
%   Optional param gather specifies whether the output vector should be
%   gathered to local lab.  
%   GATHER = 0 will not gather
%   GATHER = 1 will gather the results of forwards or adjoint
%   multiplication.
%   GATHER = 2 will gather only in forward mode.
%   GATHER = 3 will gather only in backwards (adjoint) mode.
%
%   ex.
%       M = oppKron2Lo(opDFT(15),rand(20,30));
%       x = distributed.rand(30,15);
%       y = M*x(:);
%
%   note that the second operator is applied to x and then the first
%   operator to the transpose of the result, and so x should be of
%   dimemsions [cols(op2),cols(op1)], and vectorized after distribution so
%   it is distributed along the columns evenly.
%
%   *Now oppKron2Lo also supports local x vectors and distributes them
%   before calculation( this distribution will be faster than using matlabs
%   (:) function).
%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        tflag = 0;
        gather = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron2Lo(varargin)
            op= op@opKron(varargin(1:2));
            if nargin > 2
                op.gather = varargin{end};
            end
            op.sweepflag = true;
        end % Constructor
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str=['pKron(',char(op.children{1})];
            
            % Get operators
            for i=2:length(op.children)
                str=strcat(str,[', ',char(op.children{i})]);
            end
            str=strcat(str,')');
            if op.tflag	 	
                str = strcat(str, '''');
            end
        end % Char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes is overloaded so as to call multiplication on a
        % distributed array. This multiplication will do the expected 2D
        % transform on 'x'.
        % For the moment mtimes is only implemented for right
        % multiplication
        function y=mtimes(op,x)
            if ~isa(op,'oppKron2Lo')
                error('Left multiplication not taken in account')
            elseif isnumeric(x) || isdistributed(x)
                assert( isvector(x) , 'Please use vectorized matrix')
                op.counter.plus1(op.tflag + 1 );
                y=op.multiply(x, 1 ); %use tflag to determine mode within multiply
            elseif isa(x,'opSpot')
                y = opFoG(op,x);
            else
                error(['unsupported data type: ' class(x)]);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose is overloaded to avoid wrapping the operator in an
        % opTranspose.
        function y = transpose(op)
            [m,n] = size(op);
            y = op;
            y.m = n;
            y.n = m;
            y.tflag =  ~op.tflag;
            y.permutation = op.permutation(end:-1:1);
        end
        function y = ctranspose(op)
            y = transpose(op);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double is overloaded to use a vectorized identity matrix
        function y = double(op, enable)
            if nargin < 2 || ~enable
                error(['oppKron2Lo is intended for large applications.' ...
                    ' The explicit representation will likely be very large,'...
                    ' use     double(x,1)  to proceed anyway']);
            end
            y = double(kron(op.children{1},op.children{2}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn/rrandn/zeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn is overloaded to create a distributed random vector
        function y = drandn(op)
            if ~op.tflag
                dims = [op.children{2}.n, op.children{1}.n];
            else
                dims = [op.children{2}.m, op.children{1}.m];
            end
            y = distrandnvec( dims );
        end
        function y = rrandn(op)
            if ~op.tflag
                dims = [op.children{2}.m, op.children{1}.m];
            else
                dims = [op.children{2}.n, op.children{1}.n];
            end
            y = distrandnvec( dims );
        end
        function y = dzeros(op)
            if ~op.tflag
                dims = [op.children{2}.n, op.children{1}.n];
            else
                dims = [op.children{2}.m, op.children{1}.m];
            end
            y = distzeros( dims );
        end
        function y = rzeros(op)
            if ~op.tflag
                dims = [op.children{2}.m, op.children{1}.m];
            else
                dims = [op.children{2}.n, op.children{1}.n];
            end
            y = distzeros( dims );
        end
        
    end% Methods
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            %The Kronecker product (KP) is applied to the right-hand matrix
            %taking in account the best order to apply the operators A and
            %B.
            
            %Operators
            A=op.children{1};
            B=op.children{2};
            
            %we could have been called through opSpot.applyMultiply, then
            %we need to pay attention to mode.
            if mode == 2 && ~op.tflag || mode ==1 && op.tflag
                A = A';
                B = B';
            end
                       
            %Size of the operators
            [rA,cA]=size(A);
            [rB,cB]=size(B);
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            perm = op.permutation; %Permutation to take in account.
            
            if perm(1)==2 %Classic multiplication order
                spmd
                    if iscodistributed(x)
                        y=getLocalPart(x);
                        local_width=length(y)/cB;
                        assert( mod(local_width,1) == 0, ...
                            ' x must be distributed along columns before vec')
                        y = reshape(y,cB,local_width);%reshape to local matrices
                    else
                        y = reshape(x,cB,cA);
                        y = codistributed(y);
                        y = getLocalPart(y);
                        local_width = size(y,2);
                    end
                    
                    partition = codistributed.build(local_width, ...
                        codistributor1d(2,codistributor1d.unsetPartition,...
                        [1,numlabs]));%partition of columns
                    
                    y=B*y;%apply B to local matrices, then transpose
                    y=y.';
                    
                    y = codistributed.build(y, codistributor1d...
                        (1,partition,[cA,rB]));
                    y = redistribute(y,codistributor1d(2));%distributed along cols
                    
                    y = getLocalPart(y);
                    y=A*y;%apply A to local matrices, then transpose
                    y=y.';
                    y=codistributed.build(y,codistributor1d(1,...
                        codistributor1d.unsetPartition,[rB,rA]));
                    y=redistribute(y,codistributor1d(2));
                    
                    %now vectorize y
                    y = getLocalPart(y);
                    local_size = numel(y);
                    partition = codistributed.build(local_size, ...
                        codistributor1d(2,codistributor1d.unsetPartition,...
                        [1,numlabs]));
                    y = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        partition, [rA*rB,1]));
                end
            else  %Inverted multiplication order 
                spmd
                    if iscodistributed(x)
                        y=getLocalPart(x);
                        local_width=length(y)/cB;
                        assert( mod(local_width,1) == 0, ...
                            'x must be distributed along columns before vec')
                        y = reshape(y,cB,local_width);%reshape to local matrices
                        partition = codistributed.build(local_width, ...
                            codistributor1d(2,codistributor1d.unsetPartition,...
                            [1,numlabs]));%partition of columns

                        y=y.'; %transpose since we're gonna apply A first
                        y = codistributed.build(y, codistributor1d...
                            (1,partition,[cA,cB]));
                        y = redistribute(y,codistributor1d(2));
                    else
                        y = reshape(x,cB,cA);
                        y = y.';
                        y = codistributed(y);
                    end
                    
                    y = getLocalPart(y);
                    y = A*y;%apply A to local matrices, then transpose
                    y = y.';
                    
                    y = codistributed.build(y,codistributor1d(1,...
                        codistributor1d.unsetPartition,[cB,rA]));
                    y=redistribute(y,codistributor1d(2));
                    y = getLocalPart(y);
                    y = B*y;%apply B to local matrices, no need to transpose
                    
                    %now vectorize y
                    local_size = numel(y);
                    partition = codistributed.build(local_size, ...
                        codistributor1d(2,codistributor1d.unsetPartition,...
                        [1,numlabs]));
                    y = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        partition, [rA*rB,1]));
                end
            end
            % if op.gather, y = gather(y); end %#ok<PROP,CPROP>
            if mode == 2 && ~op.tflag || mode ==1 && op.tflag % this is the correct adjoint case
                if op.gather == 1 || op.gather == 3
                    y = gather(y);
                end
            else % this is the forward case
                if op.gather == 1 || op.gather == 2
                    y = gather(y);
                end
            end % gather
        end % Multiply
    end %Protected methods
end % Classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% distrandnvec helper funtion to create a vectorized distributed
% normal random matrix
function y = distrandnvec( sz )
    spmd
        y = codistributed.randn(sz);
        codist = getCodistributor(y);
        part = codist.Partition;
        part = part.*sz(1);
        y = getLocalPart(y);
        y = y(:);
        y = codistributed.build(y, codistributor1d(1,part,...
            [sz(1)*sz(2),1]));
    end
end

% distzeros helper funtion to create a vectorized distributed
% zeros matrix
function y = distzeros( sz )
    spmd
        y = codistributed.zeros(sz);
        codist = getCodistributor(y);
        part = codist.Partition;
        part = part.*sz(1);
        y = getLocalPart(y);
        y = y(:);
        y = codistributed.build(y, codistributor1d(1,part,...
            [sz(1)*sz(2),1]));
    end
end
