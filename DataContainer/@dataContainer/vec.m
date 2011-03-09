function y = vec(x)
%VEC    Implicit vectorization of the data container
%
%   vec(A) returns data container A as if it is explicitly vectorized, but
%   is actually not. The sizes and representation of A will be of a
%   vector, but the actual vectorization of the data will only happen when
%   double is called on A.
%
%   See also: unvec, double, reshape

if ~x.reallyveced
    
    y       = x;
    d       = [x.dims{:}];
    y.dims  = {prod(d) 1};
    y.veced = true;
    y.ddims = 1;
end