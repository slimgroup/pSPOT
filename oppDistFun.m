classdef oppDistFun < oppSpot
    %OPPDISTFUN     Becuz u noe dis'fun!
    %
    %   Q = oppDistFun(A,S,F), where A is a distributed 2D matrix and S a 
    %   distributed vector, both distributed over the last dimension, 
    %   and F is a function handle that takes in local parts of A and S and
    %   gives a local part of the final answer. The arguments of F has to
    %   be standardized as: y = F(a,s,x,mode), where a and s corresponds to
    %   the local parts of A and S, and x is the distributed vector that 
    %   the operator is applied on.
    %   mode = 1 defines the forward mode
    %   mode = 2 defines the adjoint mode
    %   mode = 0 will return the sizes, complexity and linearity in an
    %   array in the following format: [m n cflag linflag]
    %   m and n has to be the conceptual size of Q so that the
    %   multiplication sizes would match.
    %   cflag is the complexity of this operator.
    %   linflag is the linearity of this operator.
    %
    %   Use case to keep in mind:
    %   (Ax - s)
    %
    %     A = length m x n
    %     s = length n
    % 
    %     P <---- constructs like oppFunComposite({A,s},F)
    % 
    %     F <-- function  y = functionF(A, B, x, mode)
    %            if mode==1
    %            y = A * x - B;
    %            elseif mode==2
    %            y = A' * x - conj(B);
    %            elseif mode=='0'
    %            retrun {m, n}
    %          end
    % 
    %     A distributed over last dim (2)
    %     s distributed over last dim (1)
    % 
    % 
    %     want P to do
    % 
    %     for each A(:,k), s(k), k = 1:n distributed
    %        y(k) = F(A(:,k),s(k),x(k))
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        fun;
        A;
        S;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDistFun(A,S,F)
            
            if nargin ~= 3 % Check for number of arguments
                error('There must be 3 arguments');
            end
            if matlabpool('size') == 0 % Check for matlabpool
                error('Matlabpool is not open');
            end
            
            if ~isdistributed(S) % Check for distributions
                error('S must be distributed');
            end
            
            if ~isdistributed(A)
                error('A must be distributed');
            end
            
            if ~isa(F,'function_handle') % Check for the function-handleness
                error('F must be a function handle'); % of F
            end
            
            % Extract parameters from function
            bleh = F(0);
            m = bleh(1); n = bleh(2); cflag = bleh(3); linflag = bleh(4);
            
            if ~isposintscalar(m) || ~isposintscalar(n) % check m and n
              error('Dimensions of operator must be positive integers.');
            end
            
            % Setup sizes
%             sizeA = size(A);
%             m = m*sizeA(end);
%             n = n*sizeA(end);
            
            % Construct oppCompositeFun
            op = op@oppSpot('DistFun', m, n);
            op.A = A;
            op.S = S;
            op.fun = F;
            op.cflag = cflag;
            op.linear = linflag;
            op.sweepflag = true;
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str='???';
            
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes is overloaded so as to call multiplication on a
        % distributed array. This multiplication will do the expected 2D
        % transform on 'x'.
        % For the moment mtimes is only implemented for right
        % multiplication
        function y=mtimes(op,x)
            if ~isa(op,'oppDistFun')
                error('Left multiplication not taken in account')
            else
                %assert( isvector(x) , 'Please use vectorized matrix')
                y=mtimes@opSpot(op,x);
            end
        end
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            A = op.A;
            S = op.S;
            F = op.fun;
            spmd
                % Setup local parts
                Aloc = getLocalPart(A);
                Sloc = getLocalPart(S);
                xloc = getLocalPart(x);
                % Setup the colon indexes
                idA(1:ndims(A) - 1) = {':'};
                idS(1:ndims(S) - 1) = {':'};
                % Setup the sizes
                sizA = size(Aloc);
%                 sizS = size(Sloc);
                
                % Compensate for single slices
                if ndims(Aloc) ~= ndims(A)
                    sizA(end+1) = 1;
                end
                
%                 if ndims(Sloc) ~= ndims(S)
%                     sizS(end+1) = 1;
%                 end
                                                
                % Preallocate y and apply function
                for k = 1:sizA(end)
                    y(idS{:},k) = F(Aloc(idA{:},k),Sloc(idS{:},k),xloc(idS{:},k),mode);
                end
                %y = codistributed.build(y,getCodistributor(x),'noCommunication');
                
            end % spmd
                
        end % multiply
        
    end % Protected methods
    
end % classdef























