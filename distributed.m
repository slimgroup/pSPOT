function y = distributed( x )
%DISTRIBUTED Summary of this function goes here
%   Detailed explanation goes here
distributedPaths = which('distributed','-ALL');
builtInDistributedPath = distributedPaths{end};
builtInDistributedPath(length(builtInDistributedPath)-13:length(builtInDistributedPath)) = [];
newDistributedPath = distributedPaths{1};
newDistributedPath(length(newDistributedPath)-13:length(newDistributedPath)) = [];
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    addpath(builtInDistributedPath);
    y = builtin(distributed, x);
    addpath(newDistributedPath);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    y = x;
end
