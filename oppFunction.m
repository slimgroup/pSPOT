function op = oppFunction(varargin)
%OPPFUNCTION Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppFunction_internal(varargin{:});
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    op = opFunction(varargin{:});
end
