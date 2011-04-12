classdef piCon < iCon
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        excoddims; % Explicit codistributed dimension of data
        excodpart; % Explicit codistributed partition of data
        imcoddims; % Implicit codistributed dimension of data
        imcodpart; % Implicit codistributed partition of data
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
            cod = cod{1};
            x.excoddims = cod.Dimension;
            x.excodpart = cod.Partition;
            x.imcoddims = cod.Dimension;
            x.imcodpart = cod.Partition;
            
        end % Constructor
        
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )
        
        % randn
        x = randn(varargin)
        
        % zeros
        x = zeros(varargin)
                
        % Serial to distributed converter
        x = distributed(data)
        
    end % Static methods
    
end % Classdef