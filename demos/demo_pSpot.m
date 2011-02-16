%% DEMO
% Two reasons for using pSPOT
% 1) Speed
% 2) Can't fit the whole thing into one node

%% Example: 

% Make sure matlabpool size is 2
poolsize = 2;
if matlabpool('size') ~= 2
    poolsize = matlabpool('size');
    if matlabpool('size') == 0
        matlabpool('open','2');
    else
        matlabpool('close','force');
        matlabpool('open','2');
    end
end

n     = 100;
nsrc  = 10;
nrec  = 10;
Irec  = randi(n,nrec,1);
nfreq = 5;
% set up Helmholtz and source matices for each frequency
for k = 1:nfreq
    Hc{k} = opGaussian(n,n,1);
    Qc{k} = opGaussian(n,nsrc);
end

% collect matrices for each frequency into supermatrices
H  = oppBlockDiag(Hc{:});
Q  = double(oppStack(Qc{:}));
P  = oppBlockDiag(nfreq,opRestriction(n,Irec));
F  = oppKron2Lo(opDFT(nfreq),opKron(opDirac(nrec),opDirac(nsrc)));

%
% word of caution:
% - DO NOT USE EXPLICIT MATRICES IN OPPSTACK!!
% - MAKE SURE THAT THE SIZES OF THE INDIVIDUAL INPUTS MATCH THE OPERATORS
% - Care should be taken so that explicit matrices are not created, 
% because oppBlockDiag replicates them internally

% solve Helmholtz equation for wavefield U
U  = H\Q; 

% check (accuracy is shitty because SPOT automatically uses 
% 10 LSQR iterations to invert the matrix. Avoid this by 
% overloading divide for you favourite operators) 
Qt = H*U;
norm(Q(:)-Qt(:))/norm(Q(:))

% make data by sampling the wavefield, then transpose and vec
% to generate vector.
D  = transpose(P*U);
d  = D(:);

% now we can do a FFT 
dt = F*d;

%% Applications
% Where A1 and A2 are helmholtz matrices
% y would be sources
% Also ask Tim for the stuff he does

% A nice use case for Stack and Dictionary
% 

U1 = opGaussian(1000,10,1); % Forward wave field for Frequency 1
U2 = opGaussian(1000,10,1); % Forward wave field for Frequency 2
V1 = opGaussian(1000,10,1); % Adjoint wave field for Frequency 1
V2 = opGaussian(1000,10,1); % Adjoint wave field for Frequency 2

% Extended image volume I
U = oppStack(U1,U2);
V = oppStack(V1,V2,1); % has to be gathered;
I = U*V'; 

% I will be stored as an opFog, which saves memory space because it is not
% explicit.

% Now we can access the image volume like so
w1 = opGaussian(1000,1);
w2 = opGaussian(1000,1);
w = double(oppStack(w1,w2));
z = I*w; 

matlabpool('close');
matlabpool('open',int2str(poolsize));

%% Example on <Insert process name here>

% Create x and S
spmd
    % Create x
    x = randn(600,500,2);
    xpart = 2*ones(1,numlabs);
    xcodist = codistributor1d(3,xpart,[600 500 8]);
    x = codistributed.build(x,xcodist,'noCommunication');
        
end

% Create S
S = randn(600,1);

% Run DFT on S
S = opDFT(600)*S;

% Run DFT along d1 of x
x = oppSweep(opDFT(600))*x;

% Hadamard of S and x
































