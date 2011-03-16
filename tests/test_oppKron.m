function test_suite = test_oppKron
initTestSuite;
end

function test_oppKron_5D
%%
dims    = 5;
lim     = 5;
DIMDIST = 2;
for i=1:dims       
    m    = 1;
    n    = 5;
    A{i} = opGaussian(m,n);
end
% Build oppKron
K1 = oppKron(A{:});
K2 = oppKron2Lo( opKron(A{1:end-1}),A{end},1);

% Construct N-D array x
xsize = cellfun(@size,K1.children,'UniformOutput',0);
for i = 1:length(xsize)
    xgsize{i} = xsize{i}(2);
end
xgsize = fliplr(xgsize);
x = randn(xgsize{:});

spmd
    x = codistributed(x,codistributor1d(DIMDIST));
end
dx = dataContainer(x);
dx = ivec(dx);

fprintf('oppKron2Lo: '),tic,y2 = K2*x(:);toc
fprintf('oppKron   : '),tic,y1 = K1*dx;toc

y1 = double(vec(unDistriCon(y1)));
assertElementsAlmostEqual(y1,y2);

end