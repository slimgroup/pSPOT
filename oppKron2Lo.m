function op = oppKron2Lo( varargin )
%OPPKRON2LO Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppKron2Lo_internal(varargin{:});
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    % Remove gather
    if isscalar(varargin{end}) && ~isa(varargin{end}, 'opSpot')
        varargin(end) = [];
    end
    op = opKron(varargin{:});
end
