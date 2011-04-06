function y = redistribute(x,varargin)

% Setup variables
dim   = ndims(x);
part  = [];
if length(varargin) == 1, dim = varargin{1}; end
if length(varargin) == 2, part = [varargin{2:end}]; end
data  = x.data;
sizes = size(x);

% Check variables
assert(dim <= ndims(x), ...
    'Distribution dimension must be within the dimensionality of x')
if ~isempty(part)
assert(sum(part) == size(x,dim),...
    'Sum of partition must be equal to size of distributed dimension')
end

% Redistribute
spmd
    % Build codistributor
    cod = codistributor1d(dim,part,sizes);
    
    % Redistribute
    data = redistribute(data,cod);
    
end % spmd

% Set variables
y = x;
cod = cod{1};
y.data = data;
y.excoddims = cod.Dimension;
y.excodpart = cod.Partition;