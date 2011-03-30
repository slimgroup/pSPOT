function e = end(A,k,n)

s = size(A);
s = [s ones(1,n-length(s)+1)];
if n == 1 && k == 1
   e = prod(s);
elseif n == ndims(A) || k < n
   e = s(k);
else % k == n || n ~= ndims(A)
   e = prod(s(k:end));
end