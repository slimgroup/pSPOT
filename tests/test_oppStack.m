function test_suite = test_oppStack
%test_oppStack  Unit tests for the Stack meta operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_oppStack_builtin
%%
   mA = 10; mB = 20; n = 10;
   A = opGaussian(mA,n);
   B = opBernoulli(mB,n);
   D = oppStack(A,B);
   utest(D,1);
   D = D';
   utest(D,1);
end

function test_oppStack_prod
%%
    A = opGaussian(10,10);
    B = opGaussian(20,10);
    D = oppStack(A,B,1);
    E = opStack(A,B);
    
    x = drandn(D,2);
    
    assertElementsAlmostEqual(D*x, E*x);
    
    x = rrandn(D,2);
    x2 = gather(x);
    
    assertElementsAlmostEqual(D'*x, E'*x2);
end

function test_oppStack_weights
%%
    m1 = randi(10); m2 = randi(10); n = randi(10);
    A1 = randn(m1,n);
    A2 = randn(m2,n);
    D = oppStack([m2 m1],A1,A2,1);
    A1 = opMatrix(m2*A1);
    A2 = opMatrix(m1*A2);
    E = opStack(A1,A2);
    
    x = drandn(D,2);
    
    assertElementsAlmostEqual(D*x, E*x);
    
    x = rrandn(D,2);
    x2 = gather(x);
    
    assertElementsAlmostEqual(D'*x, E'*x2);
end
    