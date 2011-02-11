function x = drandn(A)
%DRANDN  Normally distributed pseudorandom correctly distributed vector 
%        in the operator domain.
%
%   drandn(A) returns a pseudorandom vector in the domain of A with the
%   correct distribution matching with the operator A.
%
%   See also opSpot.rrandn, randn.

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
        xpart(i) = xpart(i) + child.n;
    end
    childnum = childnum + chipart(i);
end
xgsize = [A.n 1];

n = A.n;
if isreal(A)
    spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(n,1);
        x = redistribute(x,xcodist);
    end
else
    spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(n,1) + 1i*codistributed.randn(n,1);
        x = redistribute(x,xcodist);
    end
end
