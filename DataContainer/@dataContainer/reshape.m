function y = reshape(x,varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape(A,N1,N2,...,N) reshapes data container A into the dimensions
%   defined as [N1,N2,...,N]. Note that the number of elements must be
%   conserved.
%
%   If A is vectorized using matlab default vectorizing method (:),
%   reshaping A using the original dimensions will also redistribute A back
%   into its original distribution scheme
%
%   See also: unvec, vec, double

y = x;

if x.isdist
    data = y.Data;
    rdim = x.ddims; % Actual distributed dimension
    rvec = x.reallyveced;
    
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
    
    y.dims = varargin;
    y.Data = data;
    
else
    y.Data = reshape(x.Data,varargin{:});
    y.dims = varargin;    
end