classdef oppDistFun < oppSpot
    %OPPDISTFUN     Becuz u noe dis'fun!
    %
    %   Q = oppDistFun(A1,A2,...,AN,F,GATHER), where A1,A2,...,AN are data
    %   containers and F is a function handle that takes in local parts of 
    %   A1,A2,...,AN and gives a local part of the final answer. The 
    %   arguments of F has to be standardized as: 
    %   y = F(a1,a2,...,an,x,mode), where a1,a2,...,an corresponds to
    %   the local parts of A1,A2,...,AN, and x is the distributed vector 
    %   that the operator is applied on.
    %   mode = 1 defines the forward mode
    %   mode = 2 defines the adjoint mode
    %   mode = 0 will return the sizes, complexity and linearity in an
    %   array in the following format: [m n cflag linflag]
    %   m and n has to be the conceptual size of Q so that the
    %   multiplication sizes would match.
    %   cflag is the complexity of this operator.
    %   linflag is the linearity of this operator.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint 
    %            multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
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
        AS;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDistFun(varargin)
            
            % Check for number of arguments
            assert(nargin >= 3 ,'There must be at least 3 arguments');
            
            if matlabpool('size') == 0 % Check for matlabpool
                error('Matlabpool is not open');
            end
            % Setup and extract variables
            opgather = 0;
            if isscalar(varargin{end}) && isnumeric(varargin{end})
                assert(any(varargin{end} == [0 1 2 3]),...
                    'Gather must be 0,1,2 or 3')
                opgather      = varargin{end};
                varargin(end) = [];
            end
            
            % Check for and extract function handle
            assert(isa(varargin{end},'function_handle'),...
                'F must be a function handle');
            F             = varargin{end};
            varargin(end) = [];
            
            % Check for data containers
            assert(all(cellfun(@(p) isa(p,'dataContainer'),varargin)),...
                'A1,A2,...,AN must be data containers');
                        
            % Extract parameters from function
            bleh = F(0);
            m = bleh(1); n = bleh(2); cflag = bleh(3); linflag = bleh(4);
            
            if ~isposintscalar(m) || ~isposintscalar(n) % check m and n
              error('Dimensions of operator must be positive integers.');
            end
            
            % Setup sizes
            sizA = size(varargin{1});
            m    = m*sizA(end);
            n    = n*sizA(end);
                        
            % Construct oppCompositexun
            op           = op@oppSpot('DistFun', m, n);
            op.AS        = varargin;
            op.fun       = F;
            op.cflag     = cflag;
            op.linear    = linflag;
            op.sweepflag = true;
            op.gather    = opgather;
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str=strcat(char(op.fun),' of [');
            sizA = size(op.A);
            strsiz = int2str(sizA(1));
            for i=2:length(sizA)
                strsiz = strcat(strsiz,'x',int2str(sizA(i)));
            end
            str = strcat(str,strsiz,'] A with [');
            sizS = size(op.S);
            strsiz = int2str(sizS(1));
            for i=2:length(sizS)
                strsiz = strcat(strsiz,'-by-',int2str(sizS(i)));
            end
            str = strcat(str,strsiz,'] S');
            
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
            
            if ~isdistributed(x)
                error('x must be distributed');
            end
            
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
                sizS = size(Sloc);
                
                % Compensate for single slices
                if ndims(Aloc) ~= ndims(A)
                    sizA(end+1) = 1;
                end
                
                if ndims(Sloc) ~= ndims(S)
                    sizS(end+1) = 1;
                end
                
                % Setup partition
                ypart = codistributed.zeros(1,numlabs);
                                                
                % Apply function
                i = 0;  j = 0;
                for k = 1:sizA(end) % Iterate through the slices
                    y(i+1:i+sizA(1),1) = F(Aloc(idA{:},k),Sloc(idS{:},k),...
                        xloc(j+1:j+sizA(2),1),mode);
                    i = i+ sizA(1);
                    j = j + sizA(2);
                end
                ypart(labindex) = i; % set local partition
                % Build distributed y
                y = codistributed.build(y,codistributor1d(1,...
                    ypart,[sum(ypart) 1]),'noCommunication');
                
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
            end
                
        end % multiply
        
    end % Protected methods
    
end % classdef























