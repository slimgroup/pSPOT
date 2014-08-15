function op = oppBlockDiag( varargin )
%OPPBLOCKDIAG Summary of this function goes here
%   Detailed explanation goes here


if pSPOT.utils.hasDistComp % parallel installed
    op = oppBlockDiag_internal(varargin{:});
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    % Remove gather
    if isscalar(varargin{end}) && ~isa(varargin{end}, 'opSpot')
        varargin(end) = [];
    end
    op = opBlockDiag(varargin{:});
end
