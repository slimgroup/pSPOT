function y = reshape(varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape([DIM],X,N1,N2,...,N) reshapes data container X into the
%   dimensions defined as [N1,N2,...,N]. Note that the number of elements
%   must be conserved.
%
%   DIM specifies the distribution dimension you want your reshaped array
%   to be 'glued' in. If unspecified, the original distribution dimension
%   will be used.
%
%   Always keep in mind that reshape on distributed arrays always conserve
%   the elements locally on the labs, ie. There will be no communication
%   between labs. Therefore, the local parts size after reshaping has to be
%   the same locally. This is generally not an issue if you preserve the
%   size of the distributed dimension. Or some special symmetrical
%   distribution scheme is used.
%
%   See also: unvec, vec, double

% Check and extract dim
dim = 0;
if isa(varargin{2},'piCon')
    assert( isscalar(varargin{1}), 'Distributed dimension must be positive scalar')
    dim = varargin{1};
    varargin = varargin(2:end);
end

% Check and extract x and sizes
x = varargin{1};
varargin = varargin(2:end);
assert(isa(x,'piCon'), 'X must be a parallel data container')
sizes = [varargin{:}];

% Check for the collapsibility of reshape
% Do the calculation
imdims  = [x.imdims{:}];
while(imdims(end) == 1) % Strip singleton dimensions
    imdims(end) = [];
end
redims  = [varargin{:}];
j       = 1;
collapsed_chunk = [];
for i = 1:length(imdims)
    collapsed_chunk = [collapsed_chunk imdims(i)];
    if  prod(collapsed_chunk) == redims(j)
        collapsed_dims{j} = collapsed_chunk;
        j = j + 1;
        collapsed_chunk = [];
    elseif prod(collapsed_chunk) > redims(j)
        error('Reshape dimensions must be collapsed or multiples of implicit dimension');
    end
end

% Check for number of elements
assert(numel(x.data) == prod(sizes),'Number of elements must be conserved')

% Setup variables
data = x.data;
if ~dim, dim  = x.excoddims; end

% Pick the smallest dimension
dim = min(dim,length(sizes));
if prod(sizes) == sizes(1), dim = 1; end

spmd
    % Setup local parts
    data = getLocalPart(data);
    part = codistributed.zeros(1,numlabs);
    
    % Reshape
    if ~isempty(data)
        locsizes       = num2cell(sizes);
        locsizes{dim}  = [];
        data           = reshape(data,locsizes{:});
        part(labindex) = size(data,dim);
    else
        empty_size      = sizes;
        empty_size(dim) = 0;
        data = zeros(empty_size);
    end

    % Build codistributed
    cod  = codistributor1d(dim,part,sizes);
    data = codistributed.build(data,cod,'noCommunication');
end

% Set variables
% Compensate for 1D case
if length(collapsed_dims) == 1
    collapsed_dims{end + 1} = 1;
end
y = x;
y.data      = data;
cod         = cod{1};
y.excoddims = cod.Dimension;
y.excodpart = cod.Partition;
y.imcoddims = cod.Dimension; % Old distribution is obsolete
y.imcodpart = cod.Partition; % Old distribution is obsolete
y.imdims    = collapsed_dims;
y.perm      = 1:length(sizes); % Old permutation is obsolete
y.exdims    = sizes;