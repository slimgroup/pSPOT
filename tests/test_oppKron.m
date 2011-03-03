function test_suite = test_oppKron
    initTestSuite;
end

function test_oppKron_concept
%% Test for the viability of a n-operator oppKron based on oppKron2Lo
% Setup variables
m = randi(10);
n = randi(10);
o = 3;

for i=1:3
    A{i} = opGaussian(m,n);
end

% Dry run for x distributed over last dimension
K1 = opKron(A{:});
K2 = oppKron2Lo(opKron(A{1:2}),A{3},1);

% Build x
x = randn(n,n,n);
spmd
    x2 = codistributed(x,codistributor1d(3));
end

% Multiply
y1 = K1*x(:);
y2 = K2*x(:);

% Check
assertElementsAlmostEqual(y1,y2);

% Now try distributed over second dimension
K3 = oppKron2Lo(A{1},opKron(A{2:3}),1);

% Redistribute x
spmd
    x3 = codistributed(x,codistributor1d(2));
end
% Multiply
y3 = K3*x3(:);

% Check
assertElementsAlmostEqual(y1,y3);

end

