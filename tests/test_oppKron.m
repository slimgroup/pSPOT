function test_suite = test_oppKron
initTestSuite;
end

function test_oppKron_concept
%% Test for the concept oppKron
% Generate variables
clear all
clc
dims    = 3;
lim     = 10;
DIMDIST = 2;

for i=1:dims       
    m    = 2;
    n    = 2;
    A{i} = opGaussian(m,n);
end
% Build oppKron
K = oppKron(A{:},DIMDIST);

% Construct N-D array x
xsize = cellfun(@size,K.children,'UniformOutput',0);

for i = 1:length(xsize)
    xgsize{i} = xsize{i}(2);
end

x = randn(xgsize{:});

spmd
    x = codistributed(x,codistributor1d(DIMDIST));
end





end
