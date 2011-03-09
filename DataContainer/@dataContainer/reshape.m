function y = reshape(varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape([DIM],A,N1,N2,...,N) reshapes data container A into the dimensions
%   defined as [N1,N2,...,N]. Note that the number of elements must be
%   conserved.
%
%   DIM specifies the dimension which A is distributed. Defaults to DIM = 0
%   (undistributed) if unspecified.
%
%   If A is vectorized using matlab default vectorizing method (:),
%   reshaping A using the original dimensions will also redistribute A back
%   into its original distribution scheme
%
%   See also: unvec, vec, double

% Check and extract distributed dimension and sizes
dim = 0;
if isa(varargin{2},'dataContainer')
    assert(isscalar(varargin{1}),'Distributed dimension must be positive scalar')
    dim   = varargin{1};
    x     = varargin{2};
    sizes = varargin(3:end);
elseif isa(varargin{1}, 'dataContainer')
    x     = varargin{1};
    sizes = varargin(2:end);
else
    error('X must be a data container');    
end

% Check for distributed
assert(~xor(dim,x.isdist), 'Distributed status mismatch')
% assert(dim == x.ddims, 'Distributed dimension mismatch')

% Check for number of elements
s = size(x);
t = [sizes{:}];
assert(prod(s) == prod(t),'Number of elements must be conserved')

y = x;

if x.isdist
    data = y.Data;
    rdim = x.ddims; % Actual distributed dimension
    
    spmd        
        % Setup local parts
        data = getLocalPart(data);
        part = codistributed.zeros(1,numlabs);
        
        % Reshape
        if ~isempty(data)
        locsizes       = sizes;
        locsizes{dim}  = [];
        data           = reshape(data,locsizes{:});
        part(labindex) = size(data,dim);
        end
        
        % Build codistributed
        cod  = codistributor1d(dim,part,[sizes{:}]);
        data = codistributed.build(data,cod,'noCommunication');                
    end
    
    y.dims = sizes;
    y.Data = data;
    
else
    y.Data = reshape(x.Data,varargin{:});
    y.dims = sizes;    
end