classdef oppKron < oppSpot
    %OPPKRON    Kronecker tensor product to act on a data container.
    %   we need an oppKron that can be robust about the distribution
    %   dimension.
    %   oppKron(OP1,OP2,...,OPN) Where OPi are Spot operators or
    %   numeric matrices. Note that the last operator is applied to x then
    %   the rest of the operators to the tranpose of the result, ans so x
    %   should be of dimensions [cols(OP1),cols(OP2),...,cols(OPN)], and
    %   vectorized after distribution.
    %
    %   Optional parameter gather specifies whether the output vector
    %   should be gathered to the local lab.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %
    %   See also: oppKron2Lo
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            % Extract gather
            if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [];
            end
            
            % Standard checking and setup sizes
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            m = prod(m);
            n = prod(n);
            
            % Construct operator
            op = op@oppSpot('pKron', m, n);
            op.cflag     = cflag;
            op.linear    = linear;
            op.children  = opList;
            op.sweepflag = true;
            op.gather    = gather;
            
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
            
            % Check for dataconness of x
            assert(isa(x,'dataContainer'),...
                'X must be a data container')
            assert(x.isdist,'X must be distributed');
            
            % Remove implicit vectorization
            x = univec(x);
            
            % Setup variables
            ops = op.children;
            
            % Reversing order of children for intuitive indexing
            ops{:}
            ops = fliplr(ops);
            ops{:}
            
            % Setup spmd variables
            data  = x.data;            
            perm  = 1:length(ops);
            perm  = circshift(perm,[0 -1]);
            ddims = x.codist.Dimension;
            ddimsArray = 1:length(ops);
            ddimsArray = fliplr(ddimsArray);
            fdimsArray = circshift(ddimsArray,[0 +ddims-1]);
            
            spmd                
                % Loop through the children
                for i = 1:length(ops)
                    
                    % Update fdim
                    fdim = fdimsArray(i);
                    
                    % Multiply
                    disp(i), disp(size(ops{i})), disp(size(data))
                    labBarrier;
                    data = DataContainer.spmdMultiply(ops{i},data,mode);
                    
                    % Update gsize
                    gsize = size(data);
                    gsize = circshift(gsize,[0 -1]); disp(size(data))
                    
                    % Permute
                    [data,cod] = DataContainer.spmdPermute(data,perm,fdim,gsize);
                end
            end
            
            % Setup the variables
            y = x;
            y.data   = data;
            y.dims   = gsize{1};
            y.codist = cod{1};
            setHistory(y);
            y = ivec(y);

        end % Multiply
        
    end % Methods
    
end % classdef










