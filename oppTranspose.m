function op = oppTranspose( A )
%OPPTRANSPOSE Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppTranspose_internal(A);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opTranspose(A);
end
