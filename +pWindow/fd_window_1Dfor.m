function [ A, a ] = fd_window_1Dfor( m, p, h )
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

    A=sparse(n,m);
    a(1)=w1;
    for i=1:w1
        A(i,i)=1;
    end
    for w=2:p-1
        a(w)=w2;
        for i=1:w2
            A((w - 1)*w2 - h + i,(w - 1)*w0 - h + i)=1;
        end
    end
    a(p)=r+1;
    for i=0:r
        A(n-i,m-i)=1;
    end

end
