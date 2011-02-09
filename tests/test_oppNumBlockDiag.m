function test_suite = test_oppNumBlockDiag
% test_oppNumBlockDiag  Unit tests for the oppNumBlockDiag operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_oppNumBlockDiag_buitlin
warning('off','WarnDist:Wrongdistribution');
warning('off','WarnDist:Nodistribution');
warning('off','No:Input');
u = randn(5);
v = randn(7);
A = oppNumBlockDiag(opMatrix(u),opMatrix(v));
A.utest;
end


function test_oppNumBlockDiag_normal
%%
u = randn(5);
A1 = opBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3));
A2 = oppNumBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3));

% Test for normal spot ops
x = randn(15,2);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);

% Test for adjoint multiply
y1 = A1'*x;
y2 = A2'*x;
assertEqual(norm(y1-y2),0);

% Test for empty ops
A1 = opBlockDiag(opDFT(4),opMatrix(u),[],opDCT(3));
A2 = oppNumBlockDiag(opDFT(4),opMatrix(u),[],opDCT(3));
x = randn(12,2);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);

% Test for weights
A1 = opBlockDiag([1 2 3 4],opDFT(4),opMatrix(u),opDirac(3),opDCT(3));
A2 = oppNumBlockDiag([1 2 3 4],opDFT(4),opMatrix(u),opDirac(3),opDCT(3));
x = randn(15,2);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);

% Test for gather
A1 = oppBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3),0,1);
A2 = oppNumBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3),1);
x = randn(15,2);
y1 = A1*x;
y2 = A2*x;
assertEqual(y1,y2);
end % normal

function test_oppNumBlockDiag_repeating
A1 = opBlockDiag(3,opDFT(4));
A2 = oppNumBlockDiag(3,opDFT(4));
x = randn(12,2);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);
end % repeating

function test_oppNumBlockDiag_xdist
% Test for x vector distributed correctly
u = randn(5);
A1 = opBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3));
A2 = oppNumBlockDiag(opDFT(4),opMatrix(u),opDirac(3),opDCT(3));
x = randn(15,2);
spmd
    codist = codistributor1d(1,[9 6],[15 1]);
    x2 = codistributed(x,codist);
end
y1 = A1*x;
y2 = A2*x2;
assertEqual(norm(y1-y2),0);

% Test for x vector not distributed correctly
x2 = distributed(x);
y1 = A1*x;
y2 = A2*x2;
assertEqual(norm(y1-y2),0);
end % distributed x

function test_oppNumBlockDiag_3D
% Test for negative scalar weights
B = randn(5,4,3);
A = oppNumBlockDiag(-2,B);
assertEqual(A.weights,[-2;-2;-2]);

% Test for 3D Matrix undistributed
B = randn(6,5,4);
A1 = oppBlockDiag(B);
A2 = oppNumBlockDiag(B);
x = randn(20,1);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);

% Test for 3D Matrix correctly distributed
B2 = distributed(B);
A2 = oppNumBlockDiag(B2);
y2 = A2*x;
assertEqual(norm(y1-y2),0);

% Testing mode 2 for correctly distributed 3D Matrix
% with correctly distributed x vector
B = randn(6,5,4);
B2 = distributed(B);
A1 = oppBlockDiag([1 2 3 4],B,0,1);
A2 = oppNumBlockDiag([1 2 3 4],B2,1);
x = randn(24,1);
spmd
    codist = codistributor1d(1,[12 12],[24 1]);
    x2 = codistributed(x,codist);
end
y1 = A1'*x;
y2 = A2'*x2;
assertElementsAlmostEqual(y1,y2);

% Testing for 3D Matrix distributed wrongly
B = randn(6,5,4);
spmd
    codist = codistributor1d(2,[3 2],[6 5 4]);
    B2 = codistributed(B,codist);
end
A1 = oppBlockDiag(B);
A2 = oppNumBlockDiag(B2);
x = randn(20,1);
y1 = A1*x;
y2 = A2*x;
assertEqual(norm(y1-y2),0);
end % 3D

function test_oppNumBlockDiag_xvec
% Testing multidimensional x in vector form
B = randn(6,5,4);
A2 = oppNumBlockDiag([1 2 3 4],B,3,1);
x = randn(20,3);
x2 = distributed(x);
x2 = x2(:);
y2 = A2*x2;

% Testing for scalar kronification
B = randn(1,1,5);
A = oppNumBlockDiag(B,2,1);
x = randn(5,2);
x = x(:);
y = A*x;
y2 = A'*x;
warning('on','WarnDist:Wrongdistribution');
warning('on','WarnDist:Nodistribution');
warning('on','No:Input');
end % xvec
















