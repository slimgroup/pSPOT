function test_suite = test_oppDictionary
%test_opDictionary  Unit tests for the Dictionary meta operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_oppDictionary_mixed
   
   m = 10; nA = 20; nB = 20;
   A = opGaussian(m,nA);
   B = opBernoulli(m,nB);
   D = oppDictionary(A,B);
   
   x = D.drandn;
   
   y1 = A*x(1:nA) + B*x(nA+1:end);
   y2 = D*x;
   
   assertEqual(norm(y1-y2),0);
   assertFalse(dottest(D ))
   assertFalse(dottest(D'))
   
end

function test_oppDictionary_multiple

   D = oppDictionary(opDCT(10),opDFT(10),double(opDCT(10)));
   assertFalse(dottest(D))

   D = oppDictionary(opDCT(10),opDFT(10),[],double(opDCT(10)));
   assertFalse(dottest(D))
   assertTrue( all(size(D) == [10 30]) )
   
end

function test_oppDictionary_double
   G = opGaussian(3,5);
   E = opEye(3,4);
   R = randn(3,6);
   Z = opZeros(3,1);
   D = oppDictionary(G,E,[],R,Z);
   
   assertEqual( double(D), [double(G), double(E), double(R), double(Z)] )
end

function test_oppDictionary_utest
    