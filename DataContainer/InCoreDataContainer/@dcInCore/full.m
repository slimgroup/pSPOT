function y = full(x)

y = dcInCore(full(double(x)));
y.imdims = x.imdims;