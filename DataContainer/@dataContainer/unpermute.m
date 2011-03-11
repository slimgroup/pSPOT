function y = unpermute(x)

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