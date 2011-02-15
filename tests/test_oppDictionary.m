function test_suite = test_oppDictionary
%test_oppDictionary  Unit tests for the Dictionary meta operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_oppDictionary_builtin
%%
warning('off','pSpot:NoInput');
   m = 10; nA = 20; nB = 20;
   A = opGaussian(m,nA);
   B = opBernoulli(m,nB);
   D = oppDictionary(A,B);
   D.utest;
   D = D';
   D.utest;
end    
    
function test_oppDictionary_mixed
%%
   m = 10; nA = 20; nB = 20;
   A = opGaussian(m,nA);
   B = opBernoulli(m,nB);
   D = oppDictionary(A,B,1);
   
   x = drandn(D,2);
   x2 = gather(x);
   y1 = A*x2(1:nA,:) + B*x2(nA+1:end,:);
   y2 = D*x;
   
   assertEqual(y1,y2);   
end

function test_oppDictionary_multiple
%%
   % Test non-Spot operators
   D = oppDictionary(opDCT(5),opDFT(5),double(opDCT(5)));
   assertFalse(dottest(D,2))
    
   % Test empty operators
   D = oppDictionary(opDCT(5),opDFT(5),[],double(opDCT(5)));
   assertFalse(dottest(D,2))
   assertTrue( all(size(D) == [5 15]) )
   
end

function test_oppDictionary_double
%%
   G = opGaussian(3,5);
   E = opEye(3,4);
   R = randn(3,6);
   Z = opZeros(3,1);
   D = oppDictionary(G,E,[],R,Z);
   
   assertEqual( gather(double(D)), [double(G), double(E), double(R), double(Z)] )
   warning('on','pSpot:NoInput');
end
