classdef oppSweep < oppSpot
    %OPPSWEEP   Operator that sweeps across parallel vectors.
    %
    %   B = oppSweep(OP,GATHER) allows operator OP to perform 
    %   sweeping operations across vectors distributed across the last 
    %   dimension (which it is distributed), in parallel. OP has to be a 
    %   Spot operator
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppSweep(varargin)
                        
            % Settin' up the variables
            gather = 0;
            
            % Extract gather
            if isscalar(varargin{end}) && any(varargin{end} == [0 1 2 3])
                gather = varargin{end};
                varargin(end) = [];
            end
                        
            % Check for number of operators
            assert(length(varargin) == 1, 'Only one operator is supported');
            
            % Standard checking
            [opList,m,n,cflag,linear] = stdpspotchk(varargin{1});
            
            % Construct
            op = op@oppSpot('pSweep', m, n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.sweepflag= true;
            op.gather   = gather;
            
        end % constructor
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            
            str = ['pSweep(',char(op.children{1}),')'];
            
        end % Display
        
    end % methods
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Check for sizes
            if mode == 1
                if op.n ~= size(x,1)
                    error(...
               'Matrix dimensions must agree when multiplying by %s.',...
               char(op));
                end
            else
                if op.m ~= size(x,1)
                    error(...
               'Matrix dimensions must agree when multiplying by %s.',...
               char(op));
                end
            end
            
            % Setup variables
            A = op.children{1};
            opm = op.m; opn = op.n;
            nlabs = matlabpool('size');
            size_x = size(x);
            
            % Preallocate y
            size_rest = size_x(2:end);
            if mode == 1
                size_zeros = [opm size_rest];
            else
                size_zeros = [opn size_rest];
            end
            y = distributed.zeros(size_zeros);
            
            spmd
                % Setup local parts
                local_x = getLocalPart(x);
                local_y = getLocalPart(y);
                
                % Setup final codistributor
                finpart = codistributed.zeros(1,nlabs);
                if mode == 1
                    fingsize = [opm size_rest];
                else
                    fingsize = [opn size_rest];
                end
                
                % Multiply
                if ~isempty(local_x)
                    
                    % Setup partition
                    finpart(labindex) = size(local_x,length(size_x));
                    
                    % Vec local_x
                    vec_x = local_x(:);
                    % Setup kron
                    size_loc = size(local_x);
                    if mode == 1
                        kronop = opKron(opDirac(prod(size_loc(2:end))),A);
                    else
                        kronop = opKron(opDirac(prod(size_loc(2:end))),A');
                    end
                                        
                    % Multiply
                    local_y = kronop*vec_x;
                    
                    % Reshape
                    size_y = size(local_x);
                    if mode == 1
                        size_y(1) = opm;
                    else
                        size_y(1) = opn;
                    end
                    
                    local_y = reshape(local_y,size_y);
                    
                else
                    % Setup the empty vectors for empty labs
                    size_y = size_x;
                    if mode == 1
                        size_y(1) = opm;
                    else
                        size_y(1) = opn;
                    end
                    size_y(end) = 0;
                    local_y = zeros(size_y);
                end
                
                fincodist = codistributor1d(length(size_x),finpart,fingsize);
                
                % Build codistributed y
                y = codistributed.build(local_y,fincodist,'noCommunication');
            end
            
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
        
    end % Protected methods
    
end % classdef


















