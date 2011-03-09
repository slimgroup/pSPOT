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
x = x(:);
A = dataContainer(x);
K = oppKron2Lo(opDirac(n*o),opDFT(m));

y = K*A;

%% Test data container on Spot
m = 5;
n = 4;
o = 3;

x = randn(m,n,o) + 1i*randn(m,n,o);
x = x(:);
A = dataContainer(x);
F = opDFT(m*n*o);
y = F*A

%% Test vec and unvec
m = 5;
n = 4;
o = 3;
x = randn(m,n,o) + 1i*randn(m,n,o);
spmd
    x = codistributed(x,codistributor1d(2));
end
x = dataContainer(x)
x = vec(x);
double(x);
x
x = unvec(x)

%% Test distriCon
m = 5;
n = 4;
o = 3;
x = randn(m,n,o) + 1i*randn(m,n,o);
x = dataContainer(distributed(x))
q = double(x);
spmd,q,end
x = distriCon(x,2)
q = double(x);
spmd,q,end
