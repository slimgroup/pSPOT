classdef oppKron2Lo < oppSpot
    %oppKron2Lo  kronecker tensor product to act on a distributed vector.
    %
    %   oppKron2Lo(A,B, gather)A and B are Spot operators or numeric 
    %   matrices. Optional param gather specifies whether the output vector
    %   should be gathered to local lab.
    %   GATHER = 0 will not gather
    %   GATHER = 1 will gather the results of forwards or adjoint
    %   multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backwards (adjoint) mode.
    %
    %       %% Example: Defining seperable sparsity transforms over 
    %       different axes. Here we define a sparsity transform S that 
    %       performs Wavelet analysis on the first dimension
    %       and a 2D Curvelet analysis on the second & third dimension
    %             dim=[64,32,32];
    %             C = opCurvelet(dim(2),dim(3));
    %             W = opWavelet(dim(1),1);
    %             S = oppKron2Lo(C,W,1);
    %
    %       % Make a random 3d data-array
    %       D = distributed.randn(dim(1),prod(dim(2:end)));
    %
    %       % Check to see if the analysis followed by synthesis returns 
    %       the original signal
    %       norm(D(:)-S'*S*D(:))
    %
    %   note that the second operator is applied to x and then the first
    %   operator to the transpose of the result, and so x should be of
    %   dimemsions [cols(op2),cols(op1)], and vectorized after distribution
    %   so it is distributed along the columns evenly.
    %
    %   *Now oppKron2Lo also supports local x vectors and distributes them
    %   before calculation( this distribution will be faster than using 
    %   matlabs (:) function).
    %
    %   See also: oppKron2Lo.drandn, oppKron2Lo.rrandn, oppKron2Lo.dzeros,
    %   oppKron2Lo.rzeros
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        tflag = 0;
        permutation; %Permutation vector of intergers defining the order to
        %use when the operators (children) of the Kronecker product are
        %applied to a data vector.
        A; % Child operators
        B;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron2Lo(varargin)
            
            % Setup and extract gather
            if isscalar(varargin{end}) && isnumeric(varargin{end})
                assert(any(varargin{end} == [0 1 2 3]),...
                    'Gather must be 0 1 2 or 3')
                gat = varargin{end};
                varargin(end) = [];
            else
                gat = 0;
            end
            
            % Check for number of operators
            assert(length(varargin) == 2,'Must specify only 2 operators!')
            
            % Standard pSpot Check
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            opA = opList{1};
            opB = opList{2};
            spmd
                childA = opA;
                childB = opB;
            end
            
            % Construct operator
            op = op@oppSpot('pKron', prod(m), prod(n));
            op.cflag       = cflag;
            op.linear      = linear;
            op.sweepflag   = true;
            op.children    = [];
            op.A           = childA;
            op.B           = childB;
            op.gather      = gat;
            op.permutation = [1 2];
            op.opsn        = n;
            op.opsm        = m;
            
            % Evaluate the best permutation to use when a multiplication is
            % applied
            if ~(prod(m) == 0 || prod(n) == 0)                
                if (m(2)-n(2)) / (m(2)*n(2)) < (m(1)-n(1)) / (m(1)*n(1))
                    op.permutation = [2 1];
                end
            end
            
            % Setting up implicit dimensions of output vector
            op.ms = fliplr(cellfun(@(x) size(x,1),varargin)); % Flipped
            op.ns = fliplr(cellfun(@(x) size(x,2),varargin));
            
        end % Constructor
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            
            childs = {op.A{1} op.B{1}};
            str=['pKron(',char(childs{1})];
            
            % Get operators
            for i=2:length(childs)
                str=strcat(str,[', ',char(childs{i})]);
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
            if isa(x,'SeisDataContainer')
                y = mtimes(x,op,'swap');
            else
                
                if ~isa(op,'oppKron2Lo')
                    error('Left multiplication not taken in account')
                elseif isa(x,'opSpot')    
                    y = opFoG(op,x);
                elseif ~isa(x,'oppKron2Lo')
                    assert( isvector(x) , 'Please use vectorized matrix')
                    op.counter.plus1(op.tflag + 1 );
                    y=op.multiply(x, 1 ); % use tflag to determine mode
                                          % within multiply
                else
                    error(['unsupported data type: ' class(x)]);
                end
            end % catch
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose is overloaded to avoid wrapping the operator in an
        % opTranspose.
        function y = transpose(op)
            [m,n] = size(op);
            y = op; y.m = n; y.n = m;
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
                    ' The explicit representation will likely be very '...
                    ' large, \nuse double(x,1)  to proceed anyway']);
            end
            
            y = double(kron(op.A{1},op.B{1}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn/rrandn/zeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn is overloaded to create a distributed random vector
        function y = drandn(op)
        %DRANDN Random vector in operator domain
            if ~op.tflag
                dims = [op.opsn(2), op.opsn(1)];
            else
                dims = [op.opsm(2), op.opsm(1)];
            end
            y = distrandnvec( dims );
        end
        function y = rrandn(op)
            %RRANDN Random vector in operator range
            if ~op.tflag
                dims = [op.opsm(2), op.opsm(1)];
            else
                dims = [op.opsn(2), op.opsn(1)];
            end
            y = distrandnvec( dims );
        end
        function y = dzeros(op)
            %DZEROS Zero vector in operator domain
            if ~op.tflag
                dims = [op.opsn(2), op.opsn(1)];
            else
                dims = [op.opsm(2), op.opsm(1)];
            end
            y = distzeros( dims );
        end
        function y = rzeros(op)
            %RZEROS Zero vector in operator range
            if ~op.tflag
                dims = [op.opsm(2), op.opsm(1)];
            else
                dims = [op.opsn(2), op.opsn(1)];
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
            
            % The Kronecker product (KP) is applied to the right-hand matrix
            % taking in account the best order to apply the operators A and
            % B.
            
            %Operators
            childA = op.A;
            childB = op.B;
            mtflag = mode == 2 && ~op.tflag || mode == 1 && op.tflag;
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            if op.permutation(1)==2 %Classic multiplication order
                spmd
                    % Setup operators
                    % we could have been called through opSpot.applyMultiply,
                    % then we need to pay attention to mode.
                    if mtflag
                        A = childA';
                        B = childB';
                    else
                        A = childA;
                        B = childB;
                    end
                    
                    %Size of the operators
                    [rA,cA] = size(A);
                    [rB,cB] = size(B);
                    
                    if iscodistributed(x)
                        y              = getLocalPart(x);
                        loc_width      = length(y)/cB;
                        assert( mod(loc_width,1) == 0, ...
                            'x must be distributed along cols before vec')
                        % reshape to local 
                        y              = reshape(y,cB,loc_width); 
                    else
                        y              = reshape(x,cB,cA);
                        y              = codistributed(y);
                        y              = getLocalPart(y);
                    end
                    
                    if ~B.isDirac
                        y=B*y;% apply B to local matrices
                    end
                    
                    if ~A.isDirac
                        y=y.'; % Tranpose
                        % Build y distributed across rows
                        y = codistributed.build(y, codistributor1d...
                            (1,[],[cA,rB]),'noCommunication');
                        y = redistribute(y,codistributor1d(2));%distributed
                        % along cols
                        y = getLocalPart(y);
                        y = A*y;% apply A to local matrices, then transpose
                        y = y.'; % Now y is distributed across rows again
                        % Rebuild y as distributed across rows
                        y = codistributed.build(y,codistributor1d(1,...
                            [],[rB,rA]),'noCommunication');
                        % Redistribute y across columns
                        y = redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    % now vectorize y                    
                    % Setup part of columns
                    part = codistributed.zeros(1,numlabs);
                    part(labindex) = numel(y);
                    y              = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        part, [rA*rB,1]),'noCommunication');
                end
            else  %Inverted multiplication order
                spmd
                    % Setup operators
                    % we could have been called through opSpot.applyMultiply,
                    % then we need to pay attention to mode.
                    if mtflag
                        A = childA';
                        B = childB';
                    else
                        A = childA;
                        B = childB;
                    end
                    
                    %Size of the operators
                    [rA,cA] = size(A);
                    [rB,cB] = size(B);
                                        
                    if iscodistributed(x)
                        y              = getLocalPart(x);
                        loc_width      = length(y)/cB;
                        assert( mod(loc_width,1) == 0, ...
                            'x must be distributed along cols before vec')
                        y = reshape(y,cB,loc_width); % reshape to local 
                        
                        if ~A.isDirac
                            y = y.'; % transpose since A is applied first
                            % Build y distributed across rows
                            y = codistributed.build(y, codistributor1d...
                                (1,[],[cA,cB]),'noCommunication');
                            % Redistribute y across cols
                            y = redistribute(y,codistributor1d(2));
                        end
                    else
                        y = reshape(x,cB,cA);
                        if ~A.isDirac
                            y = y.';
                            y = codistributed(y);
                        end
                    end
                    
                    if ~A.isDirac
                        y = getLocalPart(y);
                        y = A*y;% apply A to local matrices, then transpose
                        y = y.';
                        
                        % Rebuild y distributed across rows
                        y = codistributed.build(y,codistributor1d(1,...
                            [],[cB,rA]),'noCommunication');
                        % Redistribute y across cols
                        y = redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    if ~B.isDirac
                        y = B*y;% apply B to local matrices, no need to 
                    end         % transpose
                    
                    % now vectorize y
                    % Setup part of columns
                    part = codistributed.zeros(1,numlabs);
                    part(labindex) = numel(y);
                    y              = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        part, [rA*rB,1]),'noCommunication');
                end
            end
            % if op.gather, y = gather(y); end %#ok<PROP,CPROP>
            if mode == 2 && ~op.tflag || mode ==1 && op.tflag % this is the
                if op.gather == 1 || op.gather == 3           % correct
                    y = gather(y);                            % adjoint
                end                                           % case
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
        y      = codistributed.randn(sz);
        codist = getCodistributor(y);
        part   = codist.Partition;
        part   = part.*sz(1);
        y      = getLocalPart(y);
        y      = y(:);
        y      = codistributed.build(y, codistributor1d(1,part,...
                 [sz(1)*sz(2),1]));
    end
end

% distzeros helper funtion to create a vectorized distributed
% zeros matrix
function y = distzeros( sz )
    spmd
        y      = codistributed.zeros(sz);
        codist = getCodistributor(y);
        part   = codist.Partition;
        part   = part.*sz(1);
        y      = getLocalPart(y);
        y      = y(:);
        y      = codistributed.build(y, codistributor1d(1,part,...
                 [sz(1)*sz(2),1]));
    end
end

