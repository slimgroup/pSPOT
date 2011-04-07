function test_suite = test_iCon_spgl1
%% Testing iCon with spgl1
A = opGaussian(50,100);
x = sprandn(100,1,.1);
b = A*x;

tic, y1 = spgl1(A,b); toc
tic, y2 = spgl1(A,iCon(b)); toc

norm(y1-x)
norm(y2-x)


end

