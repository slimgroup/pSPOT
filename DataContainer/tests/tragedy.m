function y = tragedy(a)

y = a;
y = dataContainer(randn(100));
z = dataContainer(randn(99));
z.data = z.data + z.data;
y.dims = z.dims;
y.data = 'TRAGEDY';