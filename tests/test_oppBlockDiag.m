function test_suite = test_oppBlockDiag
%test_oppBlockDiag  Unit tests for the opBlockDiag operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_oppBlockDiag_builtin
   n = randi(100); m = randi(100);
   A = opMatrix(randn(m,m));
   B = opMatrix(randn(n,n));
   D = oppBlockDiag(A,B);
   D.utest;
end

function test_oppBlockDiag_prod
   n = randi(100); m = randi(100);
   A = opMatrix(randn(m,m));
   B = opMatrix(randn(n,n));
   D = opBlockDiag(A,B);
   x = randn(m+n,1);
   assertElementsAlmostEqual( [A*x(1:m); B*x(m+1:end)], D*x )
   assertElementsAlmostEqual( [A'*x(1:m); B'*x(m+1:end)], D'*x )
end
