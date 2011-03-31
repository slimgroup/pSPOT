function [ A ] = funWindow1DtprFor( n, p, h )
%funWindow1DtprFor forward tapered windowing sparse array for partition-of-unity algorithms
%
%   [ A ] = funWindow1DtprFor( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      A = 'inverse' operator

    [ m os ys xs ] = pSPOT.pWindow.funWindowShape1D( n, p, h );

    A=sparse(m,n);
    for i=0:ys(1,2)-1
        A(ys(1,1)+i,ys(1,1)+i)=pSPOT.pWindow.Taper1Dcr(i+1,ys(1,2),h);
    end
    for w=2:p-1
        for i=0:ys(w,2)-1
            A(ys(w,1)+i,ys(w,1)+i-2*h*(w-1))=pSPOT.pWindow.Taper1Dmm(i+1,ys(w,2),h);;
        end
    end
    for i=0:ys(p,2)-1
        A(ys(p,1)+i,ys(p,1)+i-2*h*(p-1))=pSPOT.pWindow.Taper1Dcr(ys(p,2)-i,ys(p,2),h);
    end

    szs=size(A);
    assert(szs(1)==m,'funWindow1DtprFor: failed building the array');
    assert(szs(2)==n,'funWindow1DtprFor: failed building the array');

end
