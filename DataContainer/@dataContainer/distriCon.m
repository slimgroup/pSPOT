function y = distriCon(x,varargin)
%DISTRICON  Container distribution function
%
%   distriCon(X,DIMENSION,PARTITION,GLOBAL_SIZE) distributes or
%   redistributes data container X according to the parameter
%   specifications.
%
%   DIMENSION   specifies the dimension of distribution
%   PARTITION   is the partition of distribution
%   GLOBAL_SIZE is the global size of distribution

% Check and extract arguments
dim   = length(x.dims);
part  = [];
gsize = size(x);
l = length(varargin);
if l > 0, dim   = varargin{1}; end
if l > 1, part  = varargin{2}; end
if l > 2, gsize = varargin{3}; end

% Setup variables
data   = x.Data;
disted = x.isdist;
y      = x;

% Check for global sizes
assert( all(gsize == size(x)), 'Global size does not match')

spmd
    codist = codistributor1d(dim,part,gsize);
    if disted
        data = redistribute(data,codist);
    else
        data = codistributed(data,codist);
    end        
end

% Set variables
y.Data   = data;
y.isdist = true;
y.ddims  = dim;
y.odims  = x.odims;