classdef DataContainer < handle
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        dims = [];          %original dimensions of data
        collapsed = false;  %flag indicating if data is collapsed to a matrix
        veced = false;      %flag indicating if data is collapsed to a vector
        org = {};           %cell array describing which of the original dimensions is in each dimension of the data container
        Data = []
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = DataContainer(data)
            assert( isnumeric(data), 'DataContainer can only be created with numeric data' )
            obj.Data = distributed(data);
            obj.dims = size(data);
            obj.org = num2cell( 1:ndims(data) );
        end
        
        function d = Dims( obj )
            d = obj.dims;
        end
        
        function o = Org( obj )
            o = obj.org;
        end
    end
    
    
end