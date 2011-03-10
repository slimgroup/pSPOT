function x = reshape(x,varargin)
%RESHAPE    Reshape data container object to desired shape
%
%   reshape(X,N1,N2,...,N) reshapes data container A into the dimensions
%   defined as [N1,N2,...,N]. Note that the number of elements must be
%   conserved.
%
%   See also: unvec, vec, double

% Check and extract sizes
sizes = [varargin{:}];

% Check for number of elements
assert(numel(x) == prod(sizes),'Number of elements must be conserved')

% Setup variables
data = x.data;

if x.isdist
    data = x.data;
    dim  = x.codist.Dimension;
    
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
    
    x.data = data;    
else
    x.data = reshape(data,varargin{:});    
end

% Set variables
x.dims = sizes;
setHistory(x);