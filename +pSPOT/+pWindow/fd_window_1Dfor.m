function [ A a1 a2 ] = fd_window_1Dfor( n, p, h )
%fd_window_1Dfor forward windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    [ m w0 w1 w2 r0 r d0 ] = pSPOT.pWindow.window1Dparams( n, p, h );

    A=sparse(m,n);
    a1(1)=w1;
    a2(1)=w0;
    for i=1:w1
        A(i,i)=1;
    end
    for w=2:p-1
        a1(w)=w2;
        a2(w)=w0;
        for i=1:w2
            A((w - 1)*w2 - h + i,(w - 1)*w0 - h + i)=1;
        end
    end
    a1(p)=r+1;
    a2(p)=r0;
    for i=0:r
        A(m-i,n-i)=1;
    end
    szs=size(A);
    assert(sum(a1)==m&szs(1)==m,'tpr_window_1Dfor: failed building the array');
    assert(sum(a2)==n&szs(2)==n,'tpr_window_1Dfor: failed building the array');

end
