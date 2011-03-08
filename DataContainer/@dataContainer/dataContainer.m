classdef (InferiorClasses = {?opSpot}) dataContainer < handle
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        dims = [];      % Current dimensions of data
        odims = [];     % Original dimensions of data
        veced = false;  % flag indicating if data is implicitly vectorized
        reallyveced = false; % flag indicating if data is explicitly 
                             % vectorized
        isdist = false; % If data is distributed
        ddims = [];     % Distributed dimension
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
                
                % Extract distribution dimension
                spmd
                    cod = getCodistributor(data);
                end
                cod = cod{1};
                obj.ddims = cod.Dimension;
            else
                assert( isnumeric(data),...
                    'DataContainer can only be created with numeric data')
            end
            
            if isvector(data)
                obj.reallyveced = true;
                obj.veced = true;
            end
            
            obj.Data = data;
            obj.dims = size(data);
            obj.odims = obj.dims;
            
        end
                
    end
    
    
end