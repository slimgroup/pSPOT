function x = rrandn(A)
%RRANDN  Normally distributed pseudorandom correctly distributed vector 
%        in the operator range.
%
%   rrandn(A) returns a pseudorandom vector in the range of A.
%
%   See also opSpot.drandn, randn.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

children = A.children;
if ~isdistributed(children) % Distribute children
    opchildren = distributed(children);
else
    opchildren = children;
end

spmd, chicodist = getCodistributor(opchildren); end

chicodist = chicodist{1};
chipart = chicodist.Partition;
childnum = 0;
for i=1:matlabpool('size')
    xpart(i) = 0;
    for j=childnum+1:childnum+chipart(i)
        child = A.children{j};
        xpart(i) = xpart(i) + child.m;
    end
    childnum = childnum + chipart(i);
end
xgsize = [A.m 1];


m = A.m;

if isreal(A)
   spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(m,1);
        x = redistribute(x,xcodist);
    end
else
    spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(m,1) + 1i*codistributed.randn(m,1);
        x = redistribute(x,xcodist);
    end
end
