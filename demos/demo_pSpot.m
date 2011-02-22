%% DEMO FOR PSPOT
% Demonstration for SLIM group's pSpot Toolbox
%%
% There are two main reasons for using pSPOT:
%%
% 1)Speed
%%
% 2)Distributed memory
%% 
% Note: For this demo, please use an even number of matlabpool

%% Example: 
%
if matlabpool('size') == 0 % Setup matlabpool
    matlabpool('open','2');
end

n     = 100;
nsrc  = 10;
nrec  = 10;
Irec  = randi(n,nrec,1);
nfreq = 5;

%%
% Set up Helmholtz and source matrices for each frequency
for k = 1:nfreq
    Hc{k} = opGaussian(n,n,1);
    Qc{k} = opGaussian(n,nsrc);
end

%%
% Collect matrices for each frequency into supermatrices
H  = oppBlockDiag(Hc{:});
Q  = double(oppStack(Qc{:}));
P  = oppBlockDiag(nfreq,opRestriction(n,Irec));
F  = oppKron2Lo(opDFT(nfreq),opKron(opDirac(nrec),opDirac(nsrc)));
%%
% Note: Please make sure that the size of the individual inputs match the
% operators
%
% Solve Helmholtz equation for wavefield U
U  = H\Q;
%%
% Check (there is some margin of error because SPOT automatically uses 
% 10 LSQR iterations to invert the matrix. Avoid this by 
% overloading divide for you favourite operators) 
Qt = H*U;
norm(Q(:)-Qt(:))/norm(Q(:))
%%
% make data by sampling the wavefield, then transpose and vec
% to generate vector.
D  = transpose(P*U);
d  = D(:);
%%
% now we can do a FFT 
dt = F*d;

clear all

%% A nice use case for Stack and Dictionary
% Where A1 and A2 are helmholtz matrices
% and y would be sources
m = 100;
n = 10;
U1 = opGaussian(m,n,1); % Forward wave field for Frequency 1
U2 = opGaussian(m,n,1); % Forward wave field for Frequency 2
V1 = opGaussian(m,n,1); % Adjoint wave field for Frequency 1
V2 = opGaussian(m,n,1); % Adjoint wave field for Frequency 2

% Extended image volume I
U = oppStack(U1,U2);
V = oppStack(V1,V2,1); % has to be gathered;
I = U*V'; 
%%
% I will be stored as an opFog, which saves memory space because it is not
% explicit.
% Now we can access the image volume like so
w1 = opGaussian(m,1);
w2 = opGaussian(m,1);
w = double(oppStack(w1,w2));
z = I*w; 

clear all

%% Example on Seismic's Most Vanilla 
% Create x and S and setup the dimensions
n1 = randi(100);
n2 = randi(5);
n3 = 12;
spmd
    % Create 3D Matrix x, distributed along 3rd dimension
    x = randn(n1,n2,n3/numlabs);
    xpart = (n3/numlabs)*ones(1,numlabs);
    xcodist = codistributor1d(3,xpart,[n1 n2 n3]);
    x = codistributed.build(x,xcodist,'noCommunication');
        
end

F = opDFT(n1);
I = opDirac(n2*n3);

% Create S
s = randn(n1,1);
S = opDiag(F*s);

C = oppKron2Lo(I,F'*S*F) * x(:); % One line to rule them all

%% A more realistic example (By Felix Herrmann)
dim=[64,32,32];
C = opCurvelet(dim(2),dim(3));
W = opWavelet(dim(1),1);
S = oppKron2Lo(C,W',1);
%%
% Make data volume
D = distributed.randn(dim(1),prod(dim(2:end)));
 
norm(D(:)-S'*S*D(:))

% opW = opSplineWavelet(dim(1),1,dim(1),3,5);
% nm = opW([],0); % Determine the size
% W = opFunction(nm{1},nm{2},opW);











