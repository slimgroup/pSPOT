classdef oppDistFun < oppSpot
    %OPPDISTFUN     Becuz u noe dis'fun!
    %
    %   Q = oppDistFun(A1,A2,...,AN,F,GATHER), where A1,A2,...,AN are 
    %   distributed arrays and F is a function handle that takes in local 
    %   parts of A1,A2,...,AN and gives a local part of the final answer. 
    %   The arguments of F has to be standardized as: 
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
            
            % Turn off stupid warning
            warning('off','distcomp:codistributed:InvalidNumberOfLabs');
            
            % Store all ops as codistributed arrays inside cells
            ops = cell(1,length(varargin));
            for i = 1:length(varargin)
               d = varargin{i};
               spmd, ops{i} = d; end
            end
            opss = ops{1};
            % clear('varargin');
            
            % Check for stuffs
            c = opss{1};
            lastdim  = size(c);
            lastdim  = lastdim(end);
            lastpart = getCodistributor(c);
            lastpart = lastpart.Partition;
            for i = 2:length(opss)
                % Check for the consistency of the last dimension
                sc = size(opss{i});
                assert(sc(end) == lastdim,...
                  'The last dimension must be of the same length')
                
                % Check for isdistributed
                assert(iscodistributed(opss{i}),'A must be distributed')
                
                % Check for the distributed dimension
                cc = getCodistributor(opss{i});
                assert(length(sc) == cc.Dimension,...
                    'A must be distributed along the last dimension')
                
                % Check for partition
                assert(all(cc.Partition == lastpart),...
                    'Partition of distributed dimension must be the same')                
            end
                                    
            % Extract parameters from function
            bleh = F(0);
            m = bleh(1); n = bleh(2); cflag = bleh(3); linflag = bleh(4);
            
            if ~isposintscalar(m) || ~isposintscalar(n) % check m and n
              error('Dimensions of operator must be positive integers.');
            end
            
            % Setup sizes
            sizA = size(opss{1});
            m    = m*sizA(end);
            n    = n*sizA(end);
            clear opss;
                        
            % Construct oppCompositexun
            op           = op@oppSpot('DistFun', m, n);
            op.children  = ops;
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
           str = 'bleh';
            
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
            % Setup variables
            ops = op.children;
            F   = op.fun;
            
            % Check for the distribution of x
            assert(isdistributed(x),'X must be distributed')
                                                
            spmd
                % Setup local parts
                xloc = getLocalPart(x);
                for i = 1:length(ops)
                   ops{i} = getLocalPart(ops{i}); 
                end
                
                % Setup y
                sizeA = size(ops{1});
                y     = cell(1,sizeA(end));
                % Setup x size
                bleh  = F(0);
                xsize = bleh(2);
                
                % Loop over the slices and apply F
                n = 0;
                for i=1:sizeA(end)
                    for j = 1:length(ops) % Get last-dimensional slice
                       slice{j} = pSPOT.utils.ldind(ops{j},i); 
                    end
                    
                    % Get x slice
                    xslice = xloc(1 + n : xsize + n);
                    y{i} = F(slice{:},xslice,mode);
                    n = n + xsize;
                end
                
                % Stack y together
                y = vertcat(y{:});
                
                % Build y
                ypart = codistributed.zeros(1,numlabs);
                ypart(labindex) = size(y,1);
                ygsize = [sum(gather(ypart)) 1];
                ycod = codistributor1d(1,ypart,ygsize);
                y = codistributed.build(y,ycod,'noCommunication');
                
            end % spmd
            
%             spmd
%                 % Setup local parts
%                 Aloc = getLocalPart(A);
%                 Sloc = getLocalPart(S);
%                 xloc = getLocalPart(x);
%                 % Setup the colon indexes
%                 idA(1:ndims(A) - 1) = {':'};
%                 idS(1:ndims(S) - 1) = {':'};
%                 % Setup the sizes
%                 sizA = size(Aloc);
%                 sizS = size(Sloc);
%                 
%                 % Compensate for single slices
%                 if ndims(Aloc) ~= ndims(A)
%                     sizA(end+1) = 1;
%                 end
%                 
%                 if ndims(Sloc) ~= ndims(S)
%                     sizS(end+1) = 1;
%                 end
%                 
%                 % Setup partition
%                 ypart = codistributed.zeros(1,numlabs);
%                                                 
%                 % Apply function
%                 i = 0;  j = 0;
%                 for k = 1:sizA(end) % Iterate through the slices
%                     y(i+1:i+sizA(1),1) = F(Aloc(idA{:},k),Sloc(idS{:},k),...
%                         xloc(j+1:j+sizA(2),1),mode);
%                     i = i+ sizA(1);
%                     j = j + sizA(2);
%                 end
%                 ypart(labindex) = i; % set local partition
%                 % Build distributed y
%                 y = codistributed.build(y,codistributor1d(1,...
%                     ypart,[sum(ypart) 1]),'noCommunication');
%                 
%             end % spmd
%             
%             % Gather
%             if mode == 1
%                 if op.gather == 1 || op.gather == 2
%                     y = gather(y);
%                 end
%             else % mode == 2
%                 if op.gather == 1 || op.gather == 3
%                     y = gather(y);
%                 end
%             end
                
        end % multiply
        
    end % Protected methods
    
end % classdef























