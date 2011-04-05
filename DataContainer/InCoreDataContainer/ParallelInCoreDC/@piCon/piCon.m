classdef piCon < iCon
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        excod; % Explicit codistributor of data
        imcod; % Implicit codistributor of data
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
            x.excod = cod{1};
            x.imcod = cod{1};
            
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
        
        % empty
        x = empty(varargin)
        
        % Serial to distributed converter
        x = distributed(data)
        
    end % Static methods
    
end % Classdef