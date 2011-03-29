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
            x = x@dataContainer('InCore',dims,dims);
            x.data = data;
        end
        
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = protected )
        
        % Left Multiply
        function y = lmultiply(x,op,mode)
            
        end % Multiply
        
    end % Public methods

end % classdef