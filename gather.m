function y = gather( x )
%GATHER Summary of this function goes here
%   Detailed explanation goes here
gatherPaths = which('gather','-ALL');
builtInGatherPath = gatherPaths{end};
builtInGatherPath(length(builtInGatherPath)-8:length(builtInGatherPath)) = [];
newGatherPath = gatherPaths{1};
newGatherPath(length(newGatherPath)-8:length(newGatherPath)) = [];
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    addpath(builtInGatherPath);
    y = gather(x);
    addpath(newGatherPath);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    y = x;
end
