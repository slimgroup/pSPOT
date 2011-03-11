classdef dataContainer < handle
%DATACONTAINER  The Data Container Mother Class
%
%   dataContainer(DATA) returns a data container object containing the
%   explicit data as well as important properties of the data such as the
%   dimensions (implicit and original dimensions included), whether it is
%   implicitly or explicitly vectorized, distributed, and its distribution
%   dimension.
%
%   Since dataContainer is subclassed from the handle superclass, all 
%   instances copied from the original object points back to the same 
%   object, and data is not copied. any change made to one instance of the 
%   data container object will affect all copies of the object. 
%   (ie. copied by assignment, etc.)
%
%   Implicit methods can be called on the data container object without
%   explicitly affecting the data until needed (ie. multiplication time).
%
%   methods: size, vec, unvec, double, reshape, isdistributed, mtimes
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        dims    = [];    % Current dimensions of data
        perm    = [];    % Permutation of data for the current dimensions
        veced   = false; % flag indicating if data is vectorized 
        isdist  = false; % If data is distributed
        codist  = [];    % Current codistributor
        data    = [];    % Actual data for the container
        history = [];    % History of data container
        count   = 1;     % Counter for history
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function x = dataContainer(data)
            
            % Setup history
            History.dims = {}; % History of dimensions
            History.perm = {}; % History of permutation
            History.cod  = {}; % History of codistributors
            x.history    = History;
            
            if isdistributed(data)
                assert( strcmp(classUnderlying(data),'double'),...
                'DataContainer can only be created with numeric data' )
                x.isdist = true;
                
                % Extract distribution dimension
                spmd
                    cod      = getCodistributor(data);
                end
                x.codist     = cod{1};
            else
                assert( isnumeric(data),...
                    'DataContainer can only be created with numeric data')
            end
            
            % Check vectorization
            if isvector(data), x.veced  = true; end
            
            x.data = data;
            x.dims = size(data);
            x.perm = 1:length(x.dims);
            
            % Set history
            setHistory(x);
            
        end % Constructor
                                
    end
    
    methods (Access = protected)
        
        function setHistory(x)
            c               = x.count;
            History         = x.history;
            History.dims{c} = x.dims;
            History.perm{c} = x.perm;
            History.cod{c}  = x.codist;
            x.history       = History;
            x.count         = c + 1;
            
        end % setHistory
                
    end % Protected
    
    
end