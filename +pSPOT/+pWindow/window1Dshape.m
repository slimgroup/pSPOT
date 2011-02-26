function [ m yo yd xo xd ] = window1Dshape( n, p, h )
%window1Dshape forward windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    % initialize window vecs
    yo=zeros(p,1); yd=zeros(p,1); xo=zeros(p,1); xd=zeros(p,1);

    % test for minimal window size
    t=floor(n/(p*2));
    assert(h<t,'window1Dshape: half-overlap (%d) too large for local window size (%d). Args: (%d,%d,%d)\n',h,2*t,n,p,h);

    % get window params for target (y) oriented distribution
    m=n+(p-1)*2*h;
    f=floor(m/p); c=ceil(m/p); r=mod(m,p);

    % get window sizes and origins for target (y)
    yd(1:p)=c;
    yd(r+1:p)=f;
    yo(1)=1;
    for i=1:p-1; yo(i+1)=yo(i)+yd(i); end;

    % get window sizes and origins for source (x)
    xd(1)=yd(1)-h;
    for i=2:p-1; xd(i)=yd(i)-2*h; end
    xd(p)=yd(p)-h;
    xo(1)=1;
    for i=1:p-1; xo(i+1)=xo(i)+xd(i); end;

    % debuging only
    %disp([n sum(xd) p h m sum(yd) c r f]);
    %disp([yo']); disp([yd']); disp([xo']); disp([xd']);

    % check if everything sums up
    cn=sum(xd);
    cm=sum(yd);
    assert(cm==m,'window1Dshape: extended window sum %d != m=%d',cm,m);
    assert(cn==n,'window1Dshape: window sum %d != n=%d',cn,n);

end
