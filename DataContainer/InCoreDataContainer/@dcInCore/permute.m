function y = permute(x,varargin)
%PERMUTE    Permutation for data container
%
%   permute(X,N1,N2,...,N) permutes the data container according to the
%   order of permutation [N1,N2,...,N]
%
%   See also: unpermute

% Setup variables
perm = [varargin{:}];

% Check for permutation dimensions
assert(length(perm) == length(x.perm),'Permutation dimensions mismatch')
assert(sum(perm) == sum(x.perm),'Permutation dimensions mismatch')

% Setup future variables
% Find final distributed dimension and setup global size at the same
% time

gsize  = size(x);
gisize = isize(x);
operm  = x.perm; % Original permutation
for  i = 1:length(perm) % Find new dimension of distribution
    tgsize(i)  = gsize(perm(i));
    tgisize(i) = gisize(perm(i));
    toperm(i)  = operm(perm(i));
end

y = dcInCore(permute(x.data,perm));

y.perm   = toperm;
y.exdims = tgsize;
y.imdims = tgisize;