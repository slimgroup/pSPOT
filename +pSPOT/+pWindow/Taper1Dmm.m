function [ t ] = Taper1Dmm( i, n, h )
%Taper1Dmm 1D taper using sin() - for inner windows
%   Detailed explanation goes here

    if i<=n/2
        x=i;
    else
        x=n-i+1;
    end
    b=2*h;
    if x<b+1;
        t=sin(x*(pi/2.)/(2.*h+1.));
    else
        t=1;
    end
%    fprintf('%2d %2d %2d %f\n',b,i,x,t);
end

