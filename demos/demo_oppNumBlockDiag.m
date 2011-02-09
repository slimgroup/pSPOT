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

A = oppNumBlockDiag(A1,A2);
spmd
    x = randn(100,1); % Here x is created in the workers
    % and each worker has its own slice of x
    codist = codistributor1d(1,[100 100], [200 1]);
    xdist = codistributed.build(x,codist); % x is combined
    % as one codistributed array
end

y = A*xdist;
% x = A\y; Coming soon

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
V = oppStack(V1,V2);
I = U*V'; 

% I will be stored as an oppFog, which saves memory space because it is not
% explicit.

% Now we can access the image volume like so
z = I*randn(2000,1); 






