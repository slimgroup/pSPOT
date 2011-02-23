%% DEMO FOR PSPOT
% Demonstration for SLIM group's pSpot Toolbox
%%
% There are two main reasons for using pSPOT:
%%
% 1)Improved computational efficiency from parallel multiprocessing
%%
% 2)Ability to hadle larger dataset by distributing amongst the memory space of multiple computers
%% 
% The examples below serve to illustrate these points.
% Note: For this demo, please use an even number of matlabpool

%% Example 1: Inverting a large block-diagonal linear system with multiple data-columns
%
if matlabpool('size') == 0 % Setup matlabpool
    matlabpool('open','2');
end

n   = 1000;
n_col  = 10;
n_blocks = 6;

%%
% Set up iid Gaussian matrices and data for each block using SPOT operators
for k = 1:n_blocks
    Ac{k} = opGaussian(n,n);  % n-by-n linear system
    bc{k} = opGaussian(n,n_col);  % data
end

%%
% Collect matrices for each block into a block-diagonal system using pSPOT operators
A  = oppBlockDiag(Ac{:});
b  = double(oppStack(bc{:}));
%%
% Note: Please make sure that the size of the individual inputs match the
% operators
%
% Solve the system
xt  = A\b;
%%
% Check (there is some margin of error because SPOT automatically uses 
% 10 LSQR iterations to invert the matrix. This can be made more accurate by
% overloading divide() for your own oeprators ) 
bt = A*xt;
norm(b(:)-bt(:))/norm(b(:))

%% Example 2: Convolve every trace of a seismic datacube with a fixed signal
% Create datacube x, and fixed signal s
n1 = 128; % time samples
n2 = 20; % trace-per-gather
n3 = 12; % shot gathers
spmd
    % Create a random 3d-array x that is distributed along 3rd dimension
    x = randn(n1,n2,n3/numlabs);
    xpart = (n3/numlabs)*ones(1,numlabs);
    xcodist = codistributor1d(3,xpart,[n1 n2 n3]);
    x = codistributed.build(x,xcodist,'noCommunication');
end

% Use FFT for the convolution, create the FFT operators
F = opDFT(n1);
I = opDirac(n2*n3);

% Create fixed signal s
s = randn(n1,1);
S = opDiag(F*s);

% Create the per-trace convolution operator
C = oppKron2Lo(I,F'*S*F); % Kronecker product of a convolution operator and an identity over all the traces

% Apply the convolution by a simple multiplication
x_convolved = C * x(:);

%% Example 3: Defining seperable sparsity transforms over different axes
% Here we define a sparsity transform S that performs Wavelet analysis on the first dimension
% and a 2D Curvelet analysis on the second & third dimension
dim=[64,32,32];
C = opCurvelet(dim(2),dim(3));
W = opWavelet(dim(1),1);
S = oppKron2Lo(C,W',1);
%%
% Make a random 3d data-array
D = distributed.randn(dim(1),prod(dim(2:end)));
 
% Check to see if the analysis followed by synthesis returns the original signal
norm(D(:)-S'*S*D(:))


%% Example 4: A use-case for Stack and Dictionary operators (I CAN'T TELL WHAT THIS SECTION IS SUPPOSED TO BE ABOUT. SOMEONE PLEASE CHIME IN [TIM])
m = 100;
n = 10;
U1 = opGaussian(m,n,1); % Forward wave field for Frequency 1
U2 = opGaussian(m,n,1); % Forward wave field for Frequency 2
V1 = opGaussian(m,n,1); % Adjoint wave field for Frequency 1
V2 = opGaussian(m,n,1); % Adjoint wave field for Frequency 2

% Extended image volume I
U = oppStack(U1,U2);
V = oppStack(V1,V2,1); % has to be gathered because distributed x is not 
                       % supported for oppStack forward mode.
I = U*V'; 
%%
% I will be stored as an opFog, which saves memory space because it is not
% explicit.
% Now we can access the image volume like so
w1 = opGaussian(m,1);
w2 = opGaussian(m,1);
w = double(oppStack(w1,w2));
z = I*w; 










