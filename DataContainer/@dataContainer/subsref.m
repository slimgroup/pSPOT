function varargout = subsref(x,s)
%SUBSREF   Subscripted reference.
%
%   X.<property> returns the value stored in that property: Currently the
%   two properties that the user would be interested in is:
%   X.data    - Returns the explicit data stored in the data container
%   X.history - Returns a struct containing cell arrays of histories of the
%               dimensions, permutations and codistributor of the data 
%               container.
%
%   X(a,b,..) - where a,b,.. are indices returns the explicit 
%               elements stored within the data container as if it is a
%               Matlab array. Actually this level of subreferencing is
%               absolutely transparent. So don't expect a data container to
%               come out of this.
%
%   X(:)      - Returns a vectorized X. Note that doing this operation will
%               explicitly change all references to this object, including
%               the original, the copies and whatnot.
%
%   See also: vec, unvec, subsasgn
   
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

if length(s) > 1
   result = x;
   for i=1:length(s)
      if iscell(result)
         if strcmp(s(i).type,'{}')
            result = builtin('subsref',result,s(i));
         else
            % Apply the subsref to each element
            newresult = cell(1,length(result));
            for j=1:length(result)
               newresult{j} = subsref(result{j},s(i));
            end
            result = newresult;
         end
      else
         result = subsref(result,s(i));
      end
   end
   
   if nargout > 1
      for i=2:nargout
         varargout{i} = [];
      end
   end
   varargout{1} = result;
   
   return;
end

switch s.type
   case {'.'}
      % Set properties and flags
      varargout{1} = x.(s.subs);

   case {'{}'}
      error('Cell-indexing is not supported.');
 
   case {'()'}
       if strcmp(s.subs,':')
           varargout{1} = vec(x);
       else
           varargout{1} = subsref(x.data,s);
       end

end
