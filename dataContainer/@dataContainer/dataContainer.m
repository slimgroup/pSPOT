classdef (InferiorClasses = {?opSpot}) dataContainer < handle
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        dims = [];      % original dimensions of data
        veced = false;  % flag indicating if data is implicitly vectorized
        reallyveced = false; % flag indicating if data is explicitly 
                             % vectorized
        isdist = false; % If data is distributed
        Data = [];      % Actual data for the container
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = dataContainer(data)
            if isdistributed(data)
                assert( strcmp(classUnderlying(data),'double'),...
                'DataContainer can only be created with numeric data' )
                obj.isdist = true;
            end
            if isvector(data)
                obj.reallyveced = true;
            end
            
            obj.Data = data;
            obj.dims = size(data);
            
        end
                
    end
    
    
end