function [ B b1 b2 ] = fd_window_1Dbck( m, p, h )
%fd_window_1Dbck inverse windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    w0 = ceil(m/p);
    assert(h<w0/2,'half-overlap (%d) too large for local window size (%d)\n',h,w0);
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
        B(m-i,n-i)=1;
    end

end
