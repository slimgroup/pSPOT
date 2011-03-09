classdef (InferiorClasses = {?opSpot}) dataContainer < handle
%DATACONTAINER  The Data Container Mother Class
%
%   dataContainer(DATA) returns a data container object containing the
%   explicit data as well as important properties of the data such as the
%   dimensions (implicit and original dimensions included), whether it is
%   implicitly or explicitly vectorized, distributed, and its distribution
%   dimension.
%
%   Implicit methods can be called on the data container object without
%   explicitly affecting the data until needed (ie. multiplication time).
%
%   methods: size, vec, unvec, double, reshape, isdistributed, mtimes
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        dims        = {};    % Current dimensions of data
        odims       = {};    % Original dimensions of data
        veced       = false; % flag indicating if data is implicitly vectorized
        reallyveced = false; % flag indicating if data is explicitly 
                             % vectorized
        isdist      = false; % If data is distributed
        oddims      = 0;     % Original distributed dimension
        ddims       = 0;     % Current distributed dimension
        Data        = [];    % Actual data for the container
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
                    cod    = getCodistributor(data);
                end
                cod        = cod{1};
                obj.ddims  = cod.Dimension;
                obj.oddims = obj.ddims;
            else
                assert( isnumeric(data),...
                    'DataContainer can only be created with numeric data')
            end
            
            % Check for vectors
            if isvector(data)
                obj.reallyveced = true;
                obj.veced       = true;
            end
            
            obj.Data  = data;
            obj.dims  = num2cell(size(data));
            obj.odims = obj.dims;
            
        end % Constructor
                
    end
    
    
end