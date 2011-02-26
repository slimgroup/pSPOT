function [ B b1 b2 ] = fd_window_1Dbck( n, p, h )
%fd_window_1Dbck inverse windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    [ m w0 w1 w2 r0 r d0 ] = pSPOT.pWindow.window1Dparams( n, p, h );

    B=sparse(n,m);
    b1(1)=w0;
    b2(1)=w1;
    for i=1:w0
        B(i,i)=1;
    end
    for w=2:p-1
        b1(w)=w0;
        b2(w)=w2;
        for i=1:w0
            B((w - 1)*w0 + i,(w - 1)*w2 + i)=1;
        end
    end
    b1(p)=r0;
    b2(p)=r+1;
    for i=0:r0-1
        B(n-i,m-i)=1;
    end
    szs=size(B);
    assert(sum(b1)==n&szs(1)==n,'tpr_window_1Dbck: failed building the array');
    assert(sum(b2)==m&szs(2)==m,'tpr_window_1Dbck: failed building the array');

end
