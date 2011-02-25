function [ m w0 w1 w2 r0 r d0 ] = window1Dparams( n, p, h )
%window1Dparams forward windowing sparse array for Finite Difference algorithms
%   Detailed explanation goes here

    w0 = ceil(n/p);
    assert(h<w0/2,'window1Dparams: half-overlap (%d) too large for local window size (%d)\n',h,w0);
    w1 = w0 + h;
    w2 = w0 + 2*h;
    m = n + (p - 1)*2*h;
    r0 = mod(n, w0);
    if r0 == 0; r0 = w0; end
    r = r0 + h - 1;
    d0 = w0 - r0;

end
