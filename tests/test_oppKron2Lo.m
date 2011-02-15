function test_oppKron2Lo
%% test_opKron  Unit tests for Kronecker products
   m = matlabpool('size');
   n = m+1;
   o = n+1;
   A1 = randn(m,o) + 1i*randn(m,o);
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
