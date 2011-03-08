function obj = reshape(obj,varargin)

if obj.isdist
    data = obj.Data;
    rdim = obj.ddims; % Actual distributed dimension
    rvec = obj.reallyveced;
    
        
    spmd
        % Check for signs of redistribution
        ccod = getCodistributor(data);
        cdim = ccod.Dimension;
        if rdim   ~= cdim && rvec
            redist = true;
            dim    = cdim;
        else
            redist = false;
            dim    = rdim;
        end
        
        % Setup local parts
        data = getLocalPart(data);
        part = codistributed.zeros(1,numlabs);
        
        % Reshape
        if ~isempty(data)
        locsizes       = varargin;
        locsizes{dim}  = [];
        data           = reshape(data,locsizes{:});
        part(labindex) = size(data,dim);
        end
        
        % Build codistributed
        cod  = codistributor1d(dim,part,[varargin{:}]);
        data = codistributed.build(data,cod,'noCommunication');
        
        % Redistribute if necessary
        if redist
            data = redistribute(data,codistributor1d(rdim));            
        end
        
    end
    
    obj.dims = varargin;
    obj.Data = data;
    
else
    obj.Data = reshape(obj.Data,varargin{:});
    obj.dims = varargin;    
end