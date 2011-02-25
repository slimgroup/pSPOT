function [ A a1 a2 ] = fd_window_1Dfor( m, p, h )
%fd_window_1Dfor forward windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    w0 = ceil(m/p);
    assert(h<w0/2,'fd_window_1Dfor: half-overlap (%d) too large for local window size (%d)\n',h,w0);
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
        A(n-i,m-i)=1;
    end

end
