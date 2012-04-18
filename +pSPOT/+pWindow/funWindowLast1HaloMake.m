function [Y N PART] = funWindowLast1HaloMake( x, n, p, h )
assert(isvector(x),'Fatal error: x has to be vector')
assert(isdistributed(x),'Fatal error: x has to be distributed')
assert(matlabpool('size')==p,'Fatal error: p does not match matlabpool size')
assert(h>0,'Fatal error: half-halo has to be beger than 0')
N=n+(p-1)*2*h;
mn=prod(size(x));
assert(mod(mn,n)==0,'Fatal error: N is not a valid last dimension')
m=prod(size(x))/n;
spmd
    part=codistributor1d.defaultPartition(n);
    PART=zeros(1,p);
    PART(1)=part(1)+h;
    PART(2:end-1)=part(2:end-1)+2*h;
    PART(end)=part(end)+h;

    % get local part
    mydata = getLocalPart(x);
    % test m*part(labindex) here
    assert(prod(size(mydata))==m*part(labindex),'Fatal error: wrong local part: partition is not default?')
    mydata = reshape(mydata,[m part(labindex)]);
    otherdata=zeros([m PART(labindex)]);
    if labindex~=1 &  labindex~=numlabs
        otherdata(:,1+h:end-h)=mydata;
    elseif labindex==1
        otherdata(:,1:end-h)=mydata;
    elseif labindex==numlabs
        otherdata(:,1+h:end)=mydata;
    end

    % shift right
    if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
    if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
    %fprintf('shift right: %d <%d >%d\n',labindex,labFrom,labTo);
    halo = labSendReceive(labTo, labFrom, mydata(:,end-h+1:end));
    if labindex > 1; otherdata(:,1:1+h-1)=halo; end;

    % shift left
    if labindex>1; labTo = labindex - 1; else labTo = []; end;
    if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
    %fprintf('shift left: %d <%d >%d\n',labindex,labFrom,labTo);
    halo = labSendReceive(labTo, labFrom, mydata(:,1:1+h-1));
    if labindex < numlabs; otherdata(:,end-h+1:end)=halo; end;

    % rebuild array
    codist = codistributor1d(1,m*PART,[m*N 1]);
    Y=codistributed.build(otherdata(:),codist,'noCommunication');

end
PART = PART{1};

end
