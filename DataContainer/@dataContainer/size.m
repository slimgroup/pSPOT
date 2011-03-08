function d = size(A,dim)
%size  Dimensions of a data container
%
%   D = size(A), for a data container A, returns the dimensions of A in an
%   N-elements row vectors D = [N1,n2,...,N].
%
%   M = size(A,DIM) retuns the length of the dimension specified by
%   the scalar DIM.  Note that DIM must be within the dimensional range of 
%   A.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

d = [A.dims{:}];

if nargin == 0
   error('Not enough input arguments');

elseif nargin > 2
   error('Too many input arguments');
   
elseif nargin > 1 && ~isempty('dim')
    if nargout > 1
       error('Unknown command option.');
    end
    if dim < 1 || dim > length(A.dims)
       error('Dimension argument must be within the dimensions of A');
    else
       d = d(dim);
    end
    
end