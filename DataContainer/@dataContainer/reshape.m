function obj = reshape(obj,varargin)

if obj.isdist
    data  = obj.Data;
    ddims = obj.ddims;
    spmd
        % Setup local parts
        data = getLocalPart(data);
        part = codistributed.zeros(1,numlabs);
        
        % Reshape
        if ~isempty(data)
        locsizes = varargin;
        locsizes{ddims} = [];
        data = reshape(data,locsizes{:});
        part(labindex) = size(data,ddims);
        end
        
        % Build codistributed
        cod = codistributor1d(ddims,part,[varargin{:}]);
        data = codistributed.build(data,cod,'noCommunication');
    end
    
    obj.dims = [varargin{:}];
    obj.Data = data;
    
else
    obj.Data = reshape(obj.Data,varargin{:});
    obj.dims = [varargin{:}];    
end