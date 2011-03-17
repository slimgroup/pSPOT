%% Test for distriCon
m = 5;
n = 4;
x = randn(m,n);

x = dataContainer(x)

x = distriCon(x,1)

%% Test data reshape
m = 6;
n = 4;
o = 3;
x = randn(m,n,o) + 1i*randn(m,n,o);
x = distributed(x);
spmd
    x = redistribute(x,codistributor1d(1));
end

x = dataContainer(x)
x = reshape(2,x,9,4,2)

%% Test reshape last dimension 1
x = randn(5,4,3,2,1);
spmd
    x = codistributed(x,codistributor1d(5));
end

xd = double(unDistriCon(vec(dataContainer(x))));

assertEqual(xd,gather(x(:)));

%% Test permute & unpermute
m = 5;
n = 4;
o = 3;
p = 6;

x = randn(m,n,o,p) + 1i*randn(m,n,o,p);
x = dataContainer(distributed(x))

% Permute
x = permute(x,[1 4 3 2])
x = permute(x,[4 3 2 1])

% unpermute
x = unpermute(x)
%% Test vec and unvec
m = 5;
n = 4;
o = 3;
x = randn(m,n,o);
size(x);
spmd
    x = codistributed(x,codistributor1d(2));
end
x = dataContainer(x)
x = vec(x)

x = distriCon(x,1,[30 30]);
x = distriCon(x,1,[40 20]);
x = unvec(x)

%% Test distriCon
m = 5;
n = 4;
o = 3;
x = randn(m,n,o,p) + 1i*randn(m,n,o,p);
x = dataContainer(distributed(x))
q = double(x);
spmd,q,end
x = distriCon(x,2)
q = double(x);
spmd,q,end
