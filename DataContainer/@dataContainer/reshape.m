function x = reshape(varargin)
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
if isa(varargin{2},'dataContainer')
    assert( isscalar(varargin{1}), 'Distributed dimension must be positive scalar')
    dim = varargin{1};
    varargin = varargin(2:end);
end

% Check and extract x and sizes
x = varargin{1};
varargin = varargin(2:end);
assert(isa(x,'dataContainer'), 'X must be a data container')
sizes = [varargin{:}];

% Remove implicit
univec(x);

% Check for number of elements
assert(numel(x.data) == prod(sizes),'Number of elements must be conserved')

% Setup variables
data = x.data;

if x.isdist
%     if x.codist.Dimension > length(x.codist.Cached.GlobalSize)
%         % Codistributor bugfix where last dimension is 1. =_=
%         data = x.data;
%         if ~dim, dim  = x.codist.Dimension; end
%         spmd
%             data = getLocalPart(data);
%             part = codistributed.zeros(1,numlabs);
%             if labindex == 1
%                 data = reshape(data,sizes);
%                 part(1) = 1;
%             end
%             cod = codistributor1d(dim,part,sizes);
%             data = codistributed.build(data,cod,'noCommunication');
%         end        
%     else
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
            else
                empty_size      = sizes;
                empty_size(dim) = 0;
                data = zeros(empty_size);
            end
            
            % Build codistributed
            cod  = codistributor1d(dim,part,sizes);
            data = codistributed.build(data,cod,'noCommunication');
        end
%     end
        x.data   = data;
        x.codist = cod{1};
else
    x.data = reshape(data,varargin{:});
end

% Set variables
x.perm = 1:length(sizes); % Old permutation is void
x.dims = sizes;
if isvector(x.data), x.veced = true; end
setHistory(x);