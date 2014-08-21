function varargout = matlabpool( varargin )
%MATLABPOOL Summary of this function goes here
%   Detailed explanation goes here
matlabpoolPaths = which('matlabpool','-ALL');
builtInMatlabpoolPath = matlabpoolPaths{end};
builtInMatlabpoolPath(length(builtInMatlabpoolPath)-12:length(builtInMatlabpoolPath)) = [];
newMatlabpoolPath = matlabpoolPaths{1};
newMatlabpoolPath(length(newMatlabpoolPath)-12:length(newMatlabpoolPath)) = [];
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    addpath(builtInMatlabpoolPath);
    %[varargout{1:nargout}] = builtin('matlabpool', varargin{:});
    [varargout{1:nargout}] = matlabpool(varargin{:});
    addpath(newMatlabpoolPath);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    if nargin ~= 0
        if strcmp(varargin{1}, 'size')
            varargout{1} = 0;
        end
    end
end
