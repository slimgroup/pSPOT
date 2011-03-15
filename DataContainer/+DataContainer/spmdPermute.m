function [x,cod] = spmdPermute(x,perm,fdim,gsize)
%SPMDPERMUTE Permute function within spmd
%   *This function must be called from within an spmd block
%
%   [X,COD] = spmdPermute(X,PERM,FDIM,GSIZE)
%   Permutes codistributed array X using permutation PERM, new
%   codistributed dimension FDIM and globalsize GSIZE. Returns permutated X
%   as well as the new codistributor of x, COD.
%

% Setup local parts and Permute and re-codistribute
x    = getLocalPart(x);
part = codistributed.zeros(1,numlabs);
if ~isempty(x)
    x    = permute(x,perm);
    sizX = size(x);
    if numel(sizX) == numel(gsize) - 1
        sizX(end+1) = 1;
    end
    part(labindex) = sizX(fdim);
end
cod  = codistributor1d(fdim,part,gsize);
x    = codistributed.build(x,cod,'noCommunication');