function varargout = subsref(x,s)
%SUBSREF   Subscripted reference.
   
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
       varargout{1} = subsref(x.data,s);

end
