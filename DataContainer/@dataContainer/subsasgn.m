function x = subsasgn(x,s,b)
%SUBSASGN   Subscribed assignment.
%
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

switch s.type
   case {'.'}
        % Set properties and flags
        x.(s.subs) = b;  

   case {'{}'}
      error('Cell-index access is read-only.');
 
   case {'()'}
       x.data = subsasgn(x.data,s,b);

end
