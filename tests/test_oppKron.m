function test_suite = test_oppKron
initTestSuite;
end

function test_oppKron_5D
%% 
dims    = 5;
lim     = 10;
DIMDIST = 1;

for i=1:dims       
    m    = 3;
    n    = 4;
    A{i} = opGaussian(m,n);
end
% Build oppKron
K1 = oppKron(A{:},DIMDIST,1);
K2 = oppKron2Lo(opKron(A{1:end-1}),A{end},1);

% Construct N-D array x
xsize = cellfun(@size,K1.children,'UniformOutput',0);
for i = 1:length(xsize)
    xgsize{i} = xsize{i}(2);
end

x = randn(xgsize{:});

spmd
    x = codistributed(x,codistributor1d(DIMDIST));
end
x = x(:);


y1 = K1*x;
y2 = K2*x;

assertElementsAlmostEqual(y1,y2);

end