function x = permute(varargin)
%PERMUTE    Permutation for data container
%
%   permute(X,N1,N2,...,N) permutes the data container according to the
%   order of permutation [N1,N2,...,N]
%
%   If X is distributed, the distribution dimension will always be
%   conserved relatively. ie. no redistribution or communication of
%   elements happen under the hood.
%
%   See also: unpermute

% Setup variables
x     = varargin{1};
perm  = [varargin{2:end}];

% Check for x and distributed and permutation dimensions
assert(isa(x,'dataContainer'),'X must be a data container')
assert(length(perm) == length(x.perm),'Permutation dimensions mismatch')
assert(sum(perm) == sum(x.perm),'Permutation dimensions mismatch')

% Setup future variables
% Find final distributed dimension and setup global size at the same
% time

fdim   = 0;
gsize  = x.dims;
operm  = x.perm; % Original permutation
for  i = 1:length(perm) % Find new dimension of distribution
    if perm(i) == x.codist.Dimension
        fdim = i;
    end
    tgsize(i)  = gsize(perm(i));
    toperm(i)  = operm(perm(i));
end
gsize  = tgsize;
operm  = toperm;

if x.isdist    
    % Setup variables
    data = x.data;
        
    spmd
        % Setup local parts and Permute and re-codistribute
        data = getLocalPart(data);
        data = permute(data,perm);
        part = codistributed.zeros(1,numlabs);
        part(labindex) = size(data,fdim);
        cod  = codistributor1d(fdim,part,gsize);
        data = codistributed.build(data,cod,'noCommunication');        
        
    end % spmd
    x.data   = data;
    x.codist = cod{1};
    
else % Serial
    % Permute
    x.data  = permute(x.data,perm);
        
end
x.perm  = operm;
x.dims  = gsize;
setHistory(x);