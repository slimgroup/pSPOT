function x = rrandn(A,Ncols)
%RRANDN  Normally distributed pseudorandom correctly distributed vector 
%        in the operator range.
%
%   A.rrandn returns a pseudorandom vector in the range of A with
%   the correct distribution matching with the distribution of operator A. 
%
%   rrandn(A,NCOLS) returns the same as above except with NCOLS number of
%   columns.
%
%   Note: rrandn of oppDictionary and the transpose of oppStack will not
%   have a distributed vector returned because it does not make any sense.
%
%   See also opSpot.rrandn, randn.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
ncols = 1;
if nargin == 2 % for easy multivectoring
    ncols = Ncols;
end

if isa(A, 'oppCTranspose') || isa(A, 'oppTranspose') % Transpose mode
    A = A.children{1}; % Extract actual operator
    children = A.children;
    
    if isa(A, 'oppStack') % Stack transpose does not have 
                          % distributed rrandn
        n = A.n;
        if isreal(A)
            x = randn(n,ncols);
        else
            x = randn(n,ncols) + 1i*randn(n,ncols);
        end
        return;
    end
    
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
    xgsize = [A.n ncols];
    
    n = A.n;
    if isreal(A)
        spmd
            xcodist = codistributor1d(1,xpart,xgsize);
            x = codistributed.randn(n,ncols,codistributor1d(1));
            x = redistribute(x,xcodist);
        end
    else
        spmd
            xcodist = codistributor1d(1,xpart,xgsize);
            x = codistributed.randn(n,ncols,codistributor1d(1)) +...
                1i*codistributed.randn(n,ncols,codistributor1d(1));
            x = redistribute(x,xcodist);
        end
    end
    
    return;
end

% Forward mode

children = A.children;
if isa(A, 'oppDictionary') % Dictionary does not have distributed rrandn
    m = A.m;
    if isreal(A)
        x = randn(m,ncols);
    else
        x = randn(m,ncols) + 1i*randn(m,ncols);
    end
    return;
end

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
xgsize = [A.m ncols];

m = A.m;

if isreal(A)
   spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(m,ncols,codistributor1d(1));
        x = redistribute(x,xcodist);
    end
else
    spmd
        xcodist = codistributor1d(1,xpart,xgsize);
        x = codistributed.randn(m,ncols,codistributor1d(1)) +...
            1i*codistributed.randn(m,ncols,codistributor1d(1));
        x = redistribute(x,xcodist);
    end
end
