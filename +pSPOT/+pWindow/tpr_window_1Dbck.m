function [ A a1 a2 ] = tpr_window_1Dbck( m, p, h )
%tpr_window_1Dbck inverse tapered windowing sparse array for partition-of-unity algorithms
%   Detailed explanation goes here

    w0 = ceil(m/p);
    assert(h<w0/2,'half-overlap (%d) too large for window size (%d)\n',h,w0);
    w1 = w0 + h;
    w2 = w0 + 2*h;
    n = m + (p - 1)*2*h;
    r0 = mod(m, w0);
    if r0 == 0; r0 = w0; end
    r = r0 + h - 1;

%    fprintf('\nwindowing params:\n');
%    fprintf('\tm=%d p=%d h= %d n=%d\n',m, p, h, n);
%    fprintf('\tw0=%d/r0=%d/r=%d w1=%d w2=%d\n',w0, r0,r, w1, w2);

    A=sparse(m,n);
    a1(1)=w0;
    a2(1)=w1;
    for i=1:w1
        A(i,i)=pSPOT.pWindow.taper1Dcr(i+h,w2,h);
    end
    for w=2:p-1
        a1(w)=w0;
        a2(w)=w2;
        for i=1:w2
            A((w - 1)*w0 - h + i,(w - 1)*w2 - h + i)=pSPOT.pWindow.taper1Dmm(i,w2,h);
        end
    end
    a1(p)=r0;
    a2(p)=r+1;
    for i=0:r
        A(m-i,n-i)=pSPOT.pWindow.taper1Dcr(i+h+1,w2,h);
    end

end
