function [ F ] = funWindow1DfdFor( n, p, h )
%funWindow1DfdFor forward windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    [ m ys xs ] = pSPOT.pWindow.funWindowShape1D( n, p, h );

    F=sparse(m,n);
    for w=1:p
        for i=0:ys(w,2)-1
            F(ys(w,1)+i,ys(w,1)+i-2*h*(w-1))=1;
        end
    end

    szs=size(F);
    assert(szs(1)==m,'funWindow1DfdFor: failed building the array');
    assert(szs(2)==n,'funWindow1DfdFor: failed building the array');

end
