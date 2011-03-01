function test_suite = test_oppKron2Lo
    initTestSuite;
end
    
function test_oppKron2Lo_basic
%% test_opKron  Unit tests for Kronecker products
   m = 2;
   n = 3;
   o = 4;
   A1 = randn(n,o) + 1i*randn(n,o);
   A2 = randn(n,m) + 1i*randn(n,m);
   A3 = randn(m,m) + 1i*randn(m,m);
   A  = kron(A1,kron(A2,A3));
   BB = opKron(opMatrix(A2),opMatrix(A3));
   B  = oppKron2Lo(opMatrix(A1),BB,1);
   x  = randn(size(A,1),1) + 1i*randn(size(A,1),1);
   y  = randn(size(A,2),1) + 1i*randn(size(A,2),1);
   assertElementsAlmostEqual(A *y, B *y)
   assertElementsAlmostEqual(A'*x, B'*x)
   assertElementsAlmostEqual(A ,double(B,1))
   Bleh = double(B,1);
   assertElementsAlmostEqual(A',Bleh')
end

function test_oppKron2Lo_emptylabs
%% Test for empty labs
% Setup x
spmd
    x = codistributed.randn(100,1,1);
    xpart = [1 zeros(1,numlabs-1)];
    xgsize = [100 1 1];
    xcodist = codistributor1d(3,xpart,xgsize);
    x = redistribute(x,xcodist);
end

A = opDFT(100);
size_x = size(x);
K = oppKron2Lo(opDirac(1),A);

xvec = x(:);
y = K*xvec;
end
