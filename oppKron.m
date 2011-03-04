classdef oppKron < oppSpot
    %OPPKRON    Kronecker tensor product to act on a distributed vector.
    %   we need an oppKron that can be robust about the distribution
    %   dimension.
    %   oppKron(OP1,OP2,...,OPN,DIMDIST,GATHER) Where OPi are Spot operators or
    %   numeric matrices. Note that the last operator is applied to x then
    %   the rest of the operators to the tranpose of the result, ans so x
    %   should be of dimensions [cols(OP1),cols(OP2),...,cols(OPN)], and
    %   vectorized after distribution.
    %
    %   You have to specify the dimension at which x is distributed via the
    %   DIMDIST parameter.
    %
    %   Optional parameter gather specifies whether the output vector
    %   should be gathered to the local lab.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        dimdist; % The dimension at which x is distributed
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron(varargin)
            % Check Matlabpool
            if matlabpool('size') == 0
                error('Matlabpool is not on');
            end
            
            % Settin' up the variables
            gather = 0;
            
            % Extract dimdist and gather
            if isscalar(varargin{end-1}) && ~isa(varargin{end-1},'opSpot')
                dimdist = varargin{end-1};
                gather = varargin{end};
                varargin(end-1:end) = [];
            else
                if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot')
                    dimdist = varargin{end};
                    varargin(end) = [];
                else
                    error('Distribution dimension has to be specified');
                end
            end
            
            % Standard checking and setup sizes
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            m = prod(m);
            n = prod(n);
            
            % Construct operator
            op = op@oppSpot('pKron', m, n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.sweepflag= true;
            op.dimdist  = dimdist;
            op.gather   = gather;
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = 'pKron(';
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
        
    end % methods
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Check for distribution of x
            if ~isdistributed(x)
                error('x must be distributed');
            end
            
            % For the reshaping of x back into the correct size
            % Fetch the global size of x
            ops = op.children;
            DIMDIST    = op.dimdist;
            childsize  = cellfun(@size,ops,'UniformOutput',0);
            
            for i = 1:length(childsize)
                if mode == 1
                    xgsize{i} = childsize{i}(2);
                else
                    xgsize{i} = childsize{i}(1);
                end
            end
            
            % Reversing order of children for intuitive indexing
            u = length(ops);
            for v = 1:length(ops)
                temp(v) = ops(u);
                u = u - 1;
            end
            ops = temp;
            
            spmd
                % Reshape and redistribute x into N-D array
                xloc          = getLocalPart(x);
                nd            = length(xgsize);
                xlocsize      = xgsize;
                xlocsize{end} = [];
                xloc          = reshape(xloc,xlocsize{:});
                
                % Build codistributed N-D array distributed over last dim
                xpart           = codistributed.zeros(1,numlabs);
                xpart(labindex) = size(xloc,nd);
                defcodist       = codistributor1d(nd,xpart,[xgsize{:}]);
                x               = codistributed.build(xloc,defcodist,...
                    'noCommunication');
                
                % Redistribute into desired dimension
                x    = redistribute(x,codistributor1d(DIMDIST));
                xloc = getLocalPart(x);
                
                % Setup dimensional sizes and indices
                numops    = length(ops);
                dimArray  = 1:numops; % Array dimension
                permArray = circshift(dimArray,[0 -1]); % Permutation index
                
                % Loop through the children
                for i = 1:numops
                    
                    % Multiply
                    if i == DIMDIST % if distributed dimension
                        % Rebuild and redistribute to the next dimension
                        xpart(labindex) = size(xloc,1); % Partition
                        xsize = size(xloc); % Global size
                        xsize(1) = sum(gather(xpart));
                        distcodist = codistributor1d(1,xpart,xsize);
                        x = codistributed.build(xloc,distcodist,...
                            'noCommunication');
                        x = redistribute(x,codistributor1d(2));
                        xloc = getLocalPart(x);
                    end
                    
                    % Multiply recursively
                    if mode == 1
                        xloc = spot.utils.nDimsMultiply(ops{i},xloc);
                    else
                        xloc = spot.utils.nDimsMultiply(ops{i}',xloc);
                    end
                    
                    if i == DIMDIST
                        % Redistribute back to the original distributed
                        % dimension
                        xpart(labindex) = size(xloc,2);
                        xsize      = size(xloc);
                        xsize(2) = sum(gather(xpart));
                        nextcodist = codistributor1d(2,xpart,xsize);
                        x    = codistributed.build(xloc,nextcodist,...
                             'noCommunication');
                        x    = redistribute(x,codistributor1d(1));
                        xloc = getLocalPart(x);                        
                    end
                                        
                    % Permute to the next dimension
                    xloc = permute(xloc,permArray);
                    
                end
                % Re-distribute
                xpart(labindex) = size(xloc,DIMDIST);
                xgsize          = size(xloc);
                xgsize(DIMDIST) = sum(gather(xpart));
                xcodist         = codistributor1d(DIMDIST,xpart,xgsize);
                x               = codistributed.build(xloc,xcodist,...
                                'noCommunication');
            end % spmd            
            y = x(:);
            
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
        end % Multiply
        
    end % Methods
    
end % classdef










