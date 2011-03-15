function y = unpermute(x)
%UNPERMUTE    Un-permutation for data container
%
%   unpermute(X) permutes the data container back into the original
%   permutation as specified in the property, x.perm
%
%   If X is distributed, the distribution dimension will always be
%   conserved relatively. ie. no redistribution or communication of
%   elements happen under the hood.
%
%   See also: permute

% Remove implicit
univec(x);

% Setup permutation order
operm = x.perm;

for i = 1:length(operm)
    for j = 1:length(operm)
        if operm(j) == i;
            toperm(i) = j;
        end
    end    
end

y = permute(x,toperm);