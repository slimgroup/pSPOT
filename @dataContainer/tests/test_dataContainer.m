%% Test for data container
m = 5;
n = 4;
o = 3;

x = randn(m,n,o) + 1i*randn(m,n,o);

A = dataContainer(x);

%% Test data container multiplication
m = 5;
n = 4;
o = 3;

x = randn(m,n,o) + 1i*randn(m,n,o);
x = distributed(x);

A = dataContainer(x(:));
K = oppKron2Lo(opDirac(n*o),opDFT(m));

y = K*A;