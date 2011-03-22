function x = ivec(x)
%IVEC   Implicit vectorization of data container
%
%   ivec(x) implicitly vectorizes data container x. While the actual data
%   stored in the data container is not modified in any way, the size
%   function call on the data container will return the size of x as if it
%   is really vectorized.
%   While the data container is implicitly vectorized, performing any
%   modifying operations on x (eg, distriCon, reshape, permute, vec...)
%   will immediately void the implicit vectorized status of the data
%   container and perform the operation as though x is never implicitly
%   vectorized.
%   The sole exception to this rule would be the double function, which
%   will actually perform the vec function on the data itself and return
%   vectorized x if x is set to be implicitly vectorized.
%
%   See also: vec, unvec, double
%   

if ~x.veced
    x.ivec = true;
end