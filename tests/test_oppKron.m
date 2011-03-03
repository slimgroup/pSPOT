function test_suite = test_oppKron
initTestSuite;
end

function test_oppKron_concept
%% Test for the concept oppKron
% Generate variables
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

%% Run code
xorig = x;
op    = K;
mode  = 1;
x     = x(:);
% Check for distribution of x
if ~isdistributed(x)
    error('x must be distributed');
end

% For the reshaping of x back into the correct size
% Fetch the global size of x
ops = op.children;
DIMDIST    = op.dimdist;
childsize  = cellfun(@size,ops,'UniformOutput',0);

for i = 1:length(childsize)
    if mode == 1
        xgsize{i} = childsize{i}(2);
    else
        xgsize{i} = childsize{i}(1);
    end
end

% Reversing order of children for intuitive indexing
u = length(ops);
for v = 1:length(ops)
    temp(v) = ops(u);
    u = u - 1;
end
ops = temp;
%% spmd
spmd % Permute
    % Setup dimensional sizes and indices
    numops    = length(ops);
    dimArray  = 1:numops; % Array dimension
    permArray = circshift(dimArray,[0 -1]); % Permutation index
    distArray = [dimArray dimArray]; % Distributed dimension
    distArray = circshift(distArray,[0 -DIMDIST+1]);
    % So now the distributed dimension can be indexed with
    % distArray(i)
end
%% spmd 2
spmd
    
    % Setup x
    xloc = getLocalPart(x);
    xlocsize = xgsize;
    xlocsize{DIMDIST} = [];
    xpart = codistributed.zeros(1,numlabs);
    
    % First reshaping of vector x into N-D array
    xloc = reshape(xloc,xlocsize{:});
end
%% spmd 3
spmd
    % Loop through the children
    for i = 1:numops
        
        %                     % Multiply
        %                     if mode == 1
        %                         xloc = nDimsMultiply(ops{i},xloc);
        %                     else
        %                         xloc = nDimsMultiply(ops{i}',xloc);
        %                     end
        
        % Permute to the next dimension
        xloc = permute(xloc,permArray);
        
    end
    %                 % Re-distribute
    %                 xpart(labindex) = size(xloc,DIMDIST);
    %                 xgsize = size(xloc);
    %                 xgsize(DIMDIST) = sum(xpart);
    %                 xcodist = codistributor1d(DIMDIST,xpart,xgsize);
    %                 x = codistributed.build(xloc,xcodist,'noCommunication');
    %                 y = x;
end

%% Reverse engineering of (:) and reshape
clear all
clc
x = randn(3,3,3);
spmd
    xd   = codistributed(x,codistributor1d(2));
    xdl  = getLocalPart(xd);
    xvl  = xdl(:);
end

xv = xd(:);

spmd
    xvl = codistributed.build(xvl,getCodistributor(xv));
end

xv - xvl





















end
