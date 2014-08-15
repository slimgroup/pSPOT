function op = oppCTranspose( A )
%OPPCTRANSPOSE Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppCTranspose_internal(A);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opKron(A);
end
