classdef dcInCore < dataContainer
%DCINCORE   The In-Core Data Container Class
%
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % dcInCore Constructor
        function x = dcInCore(data)
            
            % Check for data
            assert(isnumeric(data),'Input data must be numeric')
            
            % Get sizes
            dims = size(data);
            
            % Construct class
            x = x@dataContainer('InCore',dims,num2cell(dims));
            x.data = data;
        end
        
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )
        x = randn(varargin);
        x = zeros(varargin);
        x = empty(varargin);
        
    end % Static methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = protected )
        
        % Left Multiply
        function y = lmultiply(x,op,mode)
            
        end % Multiply
        
    end % Protected methods

end % classdef