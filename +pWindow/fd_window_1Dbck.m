function [ B b ] = fd_window_1Dbck( m, p, h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    w0 = ceil(m/p);
    w1 = w0 + h;
    w2 = w0 + 2*h;
    n = m + (p - 1)*2*h;
    r0 = mod(m, w0);
    if r0 == 0; r0 = w0; end
    r = r0 + h - 1;

%    fprintf('\nwindowing params:\n');
%    fprintf('\tm=%d p=%d h= %d n=%d\n',m, p, h, n);
%    fprintf('\tw0=%d/r0=%d/r=%d w1=%d w2=%d\n',w0, r0,r, w1, w2);

    B=sparse(m,n);
    b(1)=w0;
    for i=1:w0
        B(i,i)=1;
    end
    for w=2:p-1
        b(w)=w2;
        for i=1:w0
            B((w - 1)*w0 + i,(w - 1)*w2 + i)=1;
        end
    end
    b(p)=r0;
    for i=0:r0-1
        B(m-i,n-i)=1;
    end

end
