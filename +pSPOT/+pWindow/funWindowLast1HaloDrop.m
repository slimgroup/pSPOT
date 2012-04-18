function [Y N PART] = funWindowLast1HaloDrop( x, n, p, h )
assert(isvector(x),'Fatal error: x has to be vector')
assert(isdistributed(x),'Fatal error: x has to be distributed')
assert(matlabpool('size')==p,'Fatal error: p does not match matlabpool size')
assert(h>0,'Fatal error: half-halo has to be beger than 0')
N=n+(p-1)*2*h;
mN=prod(size(x));
assert(mod(mN,N)==0,'Fatal error: N is not a valid last dimension')
m=prod(size(x))/N;
spmd
    part=codistributor1d.defaultPartition(n);
    PART=zeros(1,p);
    PART(1)=part(1)+h;
    PART(2:end-1)=part(2:end-1)+2*h;
    PART(end)=part(end)+h;

    % get local part
    mydata = getLocalPart(x);
    % test m*PART(labindex) here
    assert(prod(size(mydata))==m*PART(labindex),'Fatal error: wrong local part: partition is not the default?')
    mydata = reshape(mydata,[m PART(labindex)]);
    if labindex~=1 &  labindex~=numlabs
        otherdata=mydata(:,1+h:end-h);
    elseif labindex==1
        otherdata=mydata(:,1:end-h);
    elseif labindex==numlabs
        otherdata=mydata(:,1+h:end);
    end

    % rebuild array
    codist = codistributor1d(1,m*part,[m*n 1]);
    Y=codistributed.build(otherdata(:),codist,'noCommunication');

end
PART = part{1};

end
