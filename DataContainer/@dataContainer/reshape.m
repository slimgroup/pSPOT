function x = reshape(varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape([DIM],X,N1,N2,...,N) reshapes data container A into the dimensions
%   defined as [N1,N2,...,N]. Note that the number of elements must be
%   conserved.
%   
%   DIM specifies the distribution dimension you want your reshaped array
%   to be 'glued' in. If unspecified, the original distribution dimension
%   will be used.
%
%   See also: unvec, vec, double

% Check and extract dim
dim = 0;
if isa(varargin{2},'dataContainer')
    assert( isscalar(varargin{1}), 'Distributed dimension must be positive scalar')
    dim = varargin{1};
    varargin = varargin(2:end);
end

% Check and extract x and sizes
x = varargin{1};
assert(isa(x,'dataContainer'), 'X must be a data container')
sizes = [varargin{2:end}];

% Check for number of elements
assert(numel(x.data) == prod(sizes),'Number of elements must be conserved')

% Setup variables
data = x.data;

if x.isdist
    data = x.data;
    if ~dim, dim  = x.codist.Dimension; end
    
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
        end
        
        % Build codistributed
        cod  = codistributor1d(dim,part,sizes);
        data = codistributed.build(data,cod,'noCommunication');                
    end
    
    x.data   = data;
    x.codist = cod{1};
else
    x.data = reshape(data,varargin{:});    
end

% Set variables
x.perm = 1:length(sizes); % Old permutation is void
x.dims = sizes;
setHistory(x);