function y = uplus(x)

y = piCon(uplus(double(x)));
y.imdims    = x.imdims;
y.imcoddims = x.imcoddims;
y.imcodpart = x.imcodpart;