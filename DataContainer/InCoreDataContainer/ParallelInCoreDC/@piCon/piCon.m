classdef piCon < iCon
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        codist; % Codistributor of data
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function x = piCon(data)
            
            % Check for input data type
            assert( isdistributed(data) && ...
                strcmp(classUnderlying(data),'double'),...
            ['Parallel In-Core Data Container '...
            'can only be created with distributed numeric data'] )

            % Extract distribution dimension
            spmd
                cod  = getCodistributor(data);
            end
            
            % Construct iCon
            x = x@iCon(data);
            x.codist = cod{1};
            
        end % Constructor
        
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )
        
        % randn
        x = randn(varargin);
        
        % zeros
        x = zeros(varargin);
        
        % empty
        x = empty(varargin);
        
    end % Static methods
    
end % Classdef