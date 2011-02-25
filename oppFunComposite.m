classdef oppFunComposite < oppSpot
    %OPPQ   The final missing piece in pSpot
    %
    %   Q = oppQ(S,F,m,n,cflag,linflag) where S is a composite of vectors, and F a
    %   function handle that takes in local parts of S and gives a local 
    %   part of the final answer.
    %   m and n has to be the conceptual size of Q so that the
    %   multiplication sizes would match.
    %   cflag is the complexity of this operator.
    %   linflag is the linearity of this operator.
    %   F has to take in S and x, local parts of the respective vectors.
    %
    %   For now it is assumed that the size of the final answer will be
    %   the same as the size of x.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        fun;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppFunComposite(S,F,m,n,cflag,linflag)
            
            if nargin ~= 6 % Check for number of arguments
                error('There must be 6 arguments');
            end
            if matlabpool('size') == 0 % Check for matlabpool
                error('Matlabpool is not open');
            end
            
            if ~isa(S,'Composite') % Check for the compositeness of S
                error('S must be a composite');
            end
            
            if ~isa(F,'function_handle') % Check for the function-handleness
                error('F must be a function handle'); % of F
            end
            
            % Construct oppCompositeFun
            op = op@oppSpot('FunComposite', m, n);
            op.children = {S};
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
            if ~isa(op,'oppFunComposite')
                error('Left multiplication not taken in account')
            else
                assert( isvector(x) , 'Please use vectorized matrix')
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
            if mode ~= 1
                error('Only forward mode is allowed');
            end
            S = op.children{1};
            F = op.fun;
            spmd
                % Setup local parts
                local_x = getLocalPart(x);
                                
                % Preallocate y and apply function
                y = F(S,local_x);
                
                y = codistributed.build(y,getCodistributor(x),'noCommunication');
                
            end % spmd
                
        end % multiply
        
    end % Protected methods
    
end % classdef























