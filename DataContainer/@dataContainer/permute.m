function y = permute(varargin)
%PERMUTE    Permutation for data container
%
%   permute(x,N1,N2,...,N) permutes the data container according to the
%   order of permutation [N1,N2,...,N]
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
y = x;
fdim   = 0;
gsize  = y.dims;
ogsize = y.odims; % Original dimensions
operm  = y.perm; % Original permutation
for  i = 1:length(perm)
    if perm(i) == y.ddims
        fdim = i;
    end
    tgsize{i}  = gsize{perm(i)};
    togsize{i} = ogsize{perm(i)};
    toperm(i)  = operm(perm(i));
end
gsize  = tgsize;
ogsize = togsize;
operm  = toperm;

if y.isdist    
    % Setup variables
    data = y.Data;
        
    spmd
        % Setup local parts and Permute and re-codistribute
        data = getLocalPart(data);
        data = permute(data,perm);
        part = codistributed.zeros(1,numlabs);
        part(labindex) = size(data,fdim);
        cod  = codistributor1d(fdim,part,[gsize{:}]);
        data = codistributed.build(data,cod,'noCommunication');        
        
    end % spmd
    y.Data   = data;
    y.ddims  = fdim;
    y.oddims = fdim;
    
else % Serial
    % Permute
    y.Data  = permute(y.Data,perm);
        
end
y.perm  = operm;
y.dims  = gsize;
y.odims = ogsize;