function y = attributize(x,A)

% Check sizes
assert(all(size(x) == size(A)), 'Sizes must be consistent')

y = x;
y.data = A;