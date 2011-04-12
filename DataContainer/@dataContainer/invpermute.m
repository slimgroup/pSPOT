function y = invpermute(x)
%UNPERMUTE    Un-permutation for data container
%
%   unpermute(X) permutes the data container back into the original
%   permutation as specified in the property, x.perm
%
%   See also: permute

% Check for collapsed dimensions
for i = 1:length(x.perm)
    assert(isscalar(x.perm{i}),'Cannot unpermute with collapsed dimensions');
end

% Setup permutation order
operm = [x.perm{:}];

for i = 1:length(operm)
    for j = 1:length(operm)
        if operm(j) == i;
            toperm(i) = j;
        end
    end    
end

y = permute(x,toperm);