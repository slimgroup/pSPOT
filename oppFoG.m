function op = oppFoG( A, B )
%OPPFOG Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppFoG_internal(A,B);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opFoG(A,B);
end
