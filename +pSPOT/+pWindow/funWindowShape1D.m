function [ m ys xs ] = funWindowShape1Dfor( n, p, h )
%funWindowShape1Dfor is a support function for forward oplWindow1D* operators.
%   [ M YS XS ] = funWindowShape1Dfor( N, P, H )
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      M = length of the output vector
%      YS = (p,3) vector holding start, size, end indecies
%           of the output vector in every window
%      XS = (p,3) vector holding start, size, end indecies
%           of the input vector in every window

    % initialize shape vecs
    yo=zeros(p,1); yd=zeros(p,1); ye=zeros(p,1);
    xo=zeros(p,1); xd=zeros(p,1); xe=zeros(p,1);

    % test for minimal window size
    t=floor(n/(p*2));
    assert(h<t,'funWindowShape1Dfor: half-overlap (%d) too large for local window size (%d). Args: (%d,%d,%d)\n',h,2*t,n,p,h);

    % get window params for target (y) oriented distribution
    m=n+(p-1)*2*h;
    f=floor(m/p); c=ceil(m/p); r=mod(m,p);

    % get window sizes and origins for target (y)
    yd(1:p)=c;
    yd(r+1:p)=f;
    yo(1)=1;
    for i=1:p-1; yo(i+1)=yo(i)+yd(i); end;
    ye=yo+yd-1;

    % get window sizes and origins for source (x)
    xd(1)=yd(1)-h;
    for i=2:p-1; xd(i)=yd(i)-2*h; end
    xd(p)=yd(p)-h;
    xo(1)=1;
    for i=1:p-1; xo(i+1)=xo(i)+xd(i); end;
    xe=xo+xd-1;

    % put array together
    ys = [yo yd ye];
    xs = [xo xd xe];

    % check if everything sums up
    cn=sum(xs(:,2));
    assert(cn==n,'funWindowShape1Dfor: window sum %d != n=%d',cn,n);
    cm=sum(ys(:,2));
    assert(cm==m,'funWindowShape1Dfor: extended window sum %d != m=%d',cm,m);

    % debuging only
    %disp('funWindowShape1Dfor');
    %disp([n sum(cn) p h m sum(cm) c r f]);
    %disp([ys']); disp([xs']);

end
