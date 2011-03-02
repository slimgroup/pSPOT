function [ B ] = funWindow1DfdBck( n, p, h )
%funWindow1DfdBck inverse windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    [ m ys xs ] = pSPOT.pWindow.funWindowShape1D( n, p, h );

    B=sparse(n,m);
    for w=1:p
        for i=0:xs(w,2)-1
            B(xs(w,1)+i,xs(w,1)+i+2*h*(w-1))=1;
        end
    end

    szs=size(B);
    assert(szs(1)==n,'funWindow1DfdBck: failed building the array');
    assert(szs(2)==m,'funWindow1DfdBck: failed building the array');

end
