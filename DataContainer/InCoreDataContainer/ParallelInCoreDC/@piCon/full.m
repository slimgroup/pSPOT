function y = full(x)

y = piCon(full(double(x)));
y.imdims = x.imdims;
y.imcoddims = x.imcoddims;
y.imcodpart = x.imcodpart;