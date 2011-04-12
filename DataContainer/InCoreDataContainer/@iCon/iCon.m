classdef iCon < dataContainer
%DCINCORE   The In-Core Data Container Class
%
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % iCon Constructor
        function x = iCon(data)
            
            % Check for data
            if isdistributed(data)
                assert(strcmp(classUnderlying(data),'double'),...
                    'Input data must be numeric')
            else
                assert(isnumeric(data),'Input data must be numeric')
            end
            
            % Get sizes
            dims = size(data);
            
            % Construct class
            x = x@dataContainer('InCore',dims,num2cell(dims));
            x.data    = data;
            x.perm{1} = 1:length(size(data));
        end
        
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )
        
        % randn
        x = randn(varargin);
        
        % zeros
        x = zeros(varargin);
        
    end % Static methods
        
end % classdef