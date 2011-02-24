function [ t ] = taper1D( i, n, h )
%taper1D 1D taper using sin()
%   Detailed explanation goes here

    if i<= n/2
        x=i
    else
        x=n-i+1
    end
    if x<=h+1
        t=sin((pi/2.) * x/(2.*(h+1.)))
    else
        t=1
    end
end

