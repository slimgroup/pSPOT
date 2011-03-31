function y = uminus(x)

y = dcInCore(uminus(double(x)));
y.imdims = x.imdims;