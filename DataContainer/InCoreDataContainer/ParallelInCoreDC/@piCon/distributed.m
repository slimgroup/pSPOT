function y = distributed(data)

if isa(data,'iCon')
    y = piCon(distributed(double(data)));
    y.imdims = data.imdims;
else
    y = piCon(distributed(data));
end