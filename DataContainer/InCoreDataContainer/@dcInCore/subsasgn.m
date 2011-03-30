function x = subsasgn(x,s,b)
%SUBSASGN   Subscripted assignment.
%
%   X.<property> = A sets the value stored in that property. Though allowed
%   this is highly frowned upon. You should not mess with properties.
%
%   X(a,b,..) = A where A is a Matlab array and a,b,.. are indices that 
%               sets the explicit elements stored within the data container
%               as if it is a Matlab array. Actually this level of 
%               subsassignment is absolutely transparent. So don't pass in
%               a Data Container into a subsassignment operation.
%
%   See also: vec, unvec, subsref

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

switch s.type
   case {'.'}
        % Set properties and flags
        x.(s.subs) = b;  

   case {'{}'}
      error('Cell-index access is not supported.');
 
   case {'()'}
       x.data = subsasgn(x.data,s,b);

end
