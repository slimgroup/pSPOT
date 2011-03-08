function obj = vec(obj)
%VEC    Implicit vectorization of the data container
%
%   vec(A) returns data container A as if it is explicitly vectorized, but
%   is actually not. The sizes and representation of A will be of a
%   vector, but the actual vectorization of the data will only happen when
%   double is called on A.
%
%   See also: unvec, double, reshape

if ~obj.reallyveced
    
    tmp       = obj;
    d         = [tmp.dims{:}];
    tmp.dims  = {prod(d) 1};
    tmp.veced = true;
    obj       = tmp;
    clear('tmp');
end