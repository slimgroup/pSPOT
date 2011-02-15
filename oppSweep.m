classdef oppSweep < oppSpot
    %OPPSWEEP   Operator that sweeps across parallel vectors.
    %
    %   B = oppSweep(OP,GATHER) allows operator OP to perform 
    %   sweeping operations across vectors distributed along the second 
    %   dimension, in parallel. OP has to be a Spot operator
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
            
            % Setup variables
            A = op.children{1};
            opm = op.m; opn = op.n;
            nlabs = matlabpool('size');
            ncols = size(x,2);
            
            % Preallocate y
            y = distributed.zeros(op.m,ncols);
            
            spmd
                % Setup local parts
                local_x = getLocalPart(x);
                local_y = getLocalPart(y);
                
                % Setup final codistributor
                finpart = codistributed.zeros(1,nlabs);
                if mode == 1
                    fingsize = [opm ncols];
                else
                    fingsize = [opn ncols];
                end
                
                % Multiply
                if ~isempty(local_x)
                    
                    finpart(labindex) = size(local_x,2);
                    if mode == 1
                        local_y = A*local_x;
                    else % mode 2
                        local_y = A'*local_x;
                    end
                end
                fincodist = codistributor1d(2,finpart,fingsize);
                
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


















