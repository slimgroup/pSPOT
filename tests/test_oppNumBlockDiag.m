function test_suite = test_oppNumBlockDiag
% test_oppNumBlockDiag  Unit tests for the oppNumBlockDiag operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_oppNumBlockDiag_buitlin
%%
warning('off','WarnDist:Wrongdistribution');
warning('off','WarnDist:Nodistribution');
warning('off','No:Input');
warning('off','distcomp:codistributed:norm:usingNormest');
B = randn(5,4,2);
A = oppNumBlockDiag(distributed(B));
B1 = B(:,:,1); B2 = B(:,:,2);
C{1} = opMatrix(B1); C{2} = opMatrix(B2);
A2 = oppBlockDiag(C{:});

% Hacked dottest for oppNumBlockDiag

% Initialize
tol    = sqrt(eps);
err    = 0;
kpass  = 0;
k = 5;

Ratio = [inf,0];
for i=1:k
    x   = A2.drandn;
    y   = A2.rrandn;
    z1  = (A*x)' * y;
    z2  = x' * (A'*y);
    err = max( err, abs(z1 - z2) );
    if err < tol,
        kpass = kpass + 1;
    else
        ratio = abs(z1) / abs(z2);
        if ratio < Ratio(1), Ratio(1) = ratio; end;
        if ratio > Ratio(2), Ratio(2) = ratio; end;
    end
end

if kpass < k
    error('FAILED on %d out of %d tests\n', k-kpass, k);
    fprintf('%8s maximum absolute difference of %13.9e\n','',err);
    fprintf('%8s ratio between %13.9e and %13.9e\n', ...
        '',Ratio(1),Ratio(2));
end

end

function test_oppNumBlockDiag_weights
%%
% Test for negative scalar weights
B = randn(5,4,3);
A = oppNumBlockDiag(-2,B);
assertEqual(A.weights,[-2;-2;-2]);

% Test for normal weights
B = randn(5,4,2);
A = oppNumBlockDiag([3 4],distributed(B),1);
B1 = B(:,:,1); B2 = B(:,:,2);
C{1} = opMatrix(B1); C{2} = opMatrix(B2);
A2 = oppBlockDiag([3 4],C{:},1);

x = A2.drandn;

assertElementsAlmostEqual(A*x, A2*x);

x = A2.rrandn;

assertElementsAlmostEqual(A'*x, A2'*x);

end

function test_oppNumBlockDiag_distribution
%%
% Test for 3D Matrix undistributed
B = randn(6,5,2);
C = {opMatrix(B(:,:,1)),opMatrix(B(:,:,2))};
A1 = oppBlockDiag(C{:},1);
A2 = oppNumBlockDiag(B,1);
x = A1.drandn;
y1 = A1*x;
y2 = A2*x;
assertElementsAlmostEqual(y1,y2);

% Test for 3D Matrix correctly distributed
B2 = distributed(B);
A2 = oppNumBlockDiag(B2,1);
y2 = A2*x;
assertElementsAlmostEqual(y1,y2);

end

function test_oppNumBlockDiag_xvec
%% Testing multidimensional x in vector form
B = randn(5,4,2);
A = oppNumBlockDiag(distributed(B),3,1);
B1 = B(:,:,1); B2 = B(:,:,2);
C{1} = opKron(opDirac(3), opMatrix(B1)); 
C{2} = opKron(opDirac(3), opMatrix(B2)); 
A2 = oppBlockDiag(C{:},1);

x = drandn(A2);

y = A*x;
y2 = A2*x;

assertElementsAlmostEqual(y,y2);

%% Testing for scalar kronification
B = randn(1,1,2);
A = oppNumBlockDiag(distributed(B),3,1);
B1 = B(:,:,1); B2 = B(:,:,2);
C{1} = opKron(opDirac(3), opMatrix(B1)); 
C{2} = opKron(opDirac(3), opMatrix(B2)); 
A2 = oppBlockDiag(C{:},1);

x = drandn(A2);

y = A*x;
y2 = A2*x;

assertElementsAlmostEqual(y,y2);


warning('on','WarnDist:Wrongdistribution');
warning('on','WarnDist:Nodistribution');
warning('on','No:Input');
warning('on','distcomp:codistributed:norm:usingNormest');
end % xvec















