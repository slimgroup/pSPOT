function [ A a1 a2 ] = tpr_window_1Dbck( n, p, h )
%tpr_window_1Dbck inverse tapered windowing sparse array for partition-of-unity algorithms
%   Detailed explanation goes here

    [ m w0 w1 w2 r0 r d0 ] = pSPOT.pWindow.window1Dparams( n, p, h );

    A=sparse(n,m);
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
        A(n-i,m-i)=pSPOT.pWindow.taper1Dcr(d0+i+h+1,w2,h);
    end
    szs=size(A);
    assert(sum(a1)==n&szs(1)==n,'tpr_window_1Dbck: failed building the array');
    assert(sum(a2)==m&szs(2)==m,'tpr_window_1Dbck: failed building the array');

end
