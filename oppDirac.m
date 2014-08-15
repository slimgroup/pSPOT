function op = oppDirac( A )
%OPPDIRAC Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppDirac_internal(A);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opDirac(A);
end
