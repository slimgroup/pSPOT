function [ A a1 a2 ] = tpr_window_1Dfor( n, p, h )
%tpr_window_1Dfor forward tapered windowing sparse array for partition-of-unity algorithms
%   Detailed explanation goes here

    [ m w0 w1 w2 r0 r d0 ] = pSPOT.pWindow.window1Dparams( n, p, h );

    A=sparse(m,n);
    a1(1)=w1;
    a2(1)=w0;
    for i=1:w1
        A(i,i)=pSPOT.pWindow.taper1Dcr(i+h,w2,h);
    end
    for w=2:p-1
        a1(w)=w2;
        a2(w)=w0;
        for i=1:w2
            A((w - 1)*w2 - h + i,(w - 1)*w0 - h + i)=pSPOT.pWindow.taper1Dmm(i,w2,h);
        end
    end
    a1(p)=r+1;
    a2(p)=r0;
    for i=0:r
        A(m-i,n-i)=pSPOT.pWindow.taper1Dcr(d0+i+h+1,w2,h);
    end

end
