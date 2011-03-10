function x = distriCon(x,varargin)
%DISTRICON  Container distribution function
%
%   distriCon(X,DIMENSION,PARTITION,GLOBAL_SIZE) distributes or
%   redistributes data container X according to the parameter
%   specifications.
%
%   distriCon(X,CODIST) does the same thing, except here CODIST must be a
%   codistributor1d
%
%   DIMENSION   specifies the dimension of distribution. If omitted, the
%               last dimension will be used.
%   PARTITION   is the partition of distribution. If omitted, the default
%               partition will be used.
%   GLOBAL_SIZE is the global size of distribution. If omitted, the default
%               global size will be used.

% Check and extract arguments
dim   = length(x.dims);
part  = [];
gsize = size(x);
l     = length(varargin);
if isa(varargin{1},'codistributor1d') % Codistributor case
    dim = varargin{1}.Dimension;
    part = varargin{1}.Partition;
    gsize = varargin{1}.Cached.GlobalSize;
else % Manual params
    if l > 0, dim   = varargin{1}; end
    if l > 1, part  = varargin{2}; end
    if l > 2, gsize = varargin{3}; end
end

% Setup variables
data   = x.data;
disted = x.isdist;

% Check for partition and global sizes
assert( all(gsize == size(x)), 'Global size mismatch')
if ~isempty(part)
    assert( sum(part) == gsize(dim), 'Partition does not match global size')
end

spmd
    codist   = codistributor1d(dim,part,gsize);
    if disted
        data = redistribute(data,codist);
    else
        data = codistributed(data,codist);
    end
end

% Set variables
x.data             = data;
x.isdist           = true;
x.codist           = codist{1};
x.setHistory;
