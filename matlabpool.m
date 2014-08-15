function varargout = matlabpool( varargin )
%MATLABPOOL Summary of this function goes here
%   Detailed explanation goes here
slimapps = getenv('SLIM_APPS');
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    rmpath([slimapps '/tools/utilities/pSPOT/'])
    %[varargout{1:nargout}] = builtin('matlabpool', varargin{:});
    [varargout{1:nargout}] = matlabpool(varargin{:});
    addpath([slimapps '/tools/utilities/pSPOT/'])
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    if nargin ~= 0
        if strcmp(varargin{1}, 'size')
            varargout{1} = 0;
        end
    end
end
