function y = uminus(x)

y = piCon(uminus(double(x)));
y.imdims    = x.imdims;
y.imcoddims = x.imcoddims;
y.imcodpart = x.imcodpart;