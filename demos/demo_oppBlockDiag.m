%% DEMO
% Two reasons for using pSPOT
% 1) Speed
% 2) Can't fit the whole thing into one node

%% Example: 

if matlabpool('size') == 0
    matlabpool('open','2');
end
A1 = opGaussian(100,100,1);
A2 = opGaussian(100,100,1);
% Note: Care should be taken so that explicit 
% matrices are not created, because oppNumBlockDiag
% replicates them internally

A  = oppBlockDiag(A1,A2);
x1 = opGaussian(100,1);
x2 = opGaussian(100,1);
x = double(oppStack(x1,x2));

% word of caution:
% - DO NOT USE EXPLICIT MATRICES IN OPPSTACK!!
% - MAKE SURE THAT THE SIZES OF THE INDIVIDUAL INPUTS MATCH THE OPERATORS

y  = A*x;
xt = A\y; 

% check 
norm(x-xt)/norm(x)

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






