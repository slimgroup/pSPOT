classdef dataContainer
%DATACONTAINER  The Data Container Mother Class
%
%   x = dataContainer(DATA) returns a data container object containing the
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
%   The following properties can be accessed via dot notation:
%
%   x.data returns the actual data stored in the data container. Note that
%   any changes made to this data by itself is not monitored by the data
%   container and thus may break its history.
%
%   data container methods: 
%   distriCon, reshape, vec, unvec, permute, unpermute, mtimes
%
%   overloaded Matlab methods:
%   size, norm, subsref, subsasgn
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        exdims = []; % Explicit dimensions of data
        imdims = {}; % Implicit dimensions of data
        perm   = []; % Permutation of data (since original construction)
        type   = ''; % Type of data container
        strict = false; % Strict flag for elementary operations
    end
    
    properties ( Access = protected )
        data   = []; % Actual data for the container
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % DataCon Constructor
        function x = dataContainer(type,exdims,imdims)
            
            % Check number of arguments
            assert(nargin == 3,'There must be exactly 3 input arguments')
            
            % Set attributes
            x.type   = type;
            x.exdims = exdims;
            x.imdims = imdims;
            
        end % Constructor
        
        % Access Methods
        function value = get.type(x)
            value = x.type;
        end
        
        function value = get.exdims(x)
            value = x.exdims;
        end
        
        function value = get.imdims(x)
            value = x.imdims;
        end
        
        function value = get.data(x)
            value = x.data;
        end
        
        % Set Methods
        function x = set.type(x,value)
            x.type   = value;
        end
        
        function x = set.exdims(x,value)
            x.exdims = value;
        end
        
        function x = set.imdims(x,value)
            x.imdims = value;
        end
        
        function x = set.data(x,value)
            x.data   = value;
        end
                
    end % Public methods
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pure Protected Virtual Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Abstract, Access = protected)
        
        % Left Multiply
        y = lmultiply(x,op,mode);
                
    end % Pure virtual methods
    
    
end