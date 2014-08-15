function op = oppMatrix( A, varargin )
%OPPMATRIX Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppMatrix_internal(A, varargin{:});
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opMatrix(A, varargin{:});
end
