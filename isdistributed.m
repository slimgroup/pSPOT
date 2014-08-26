function tf = isdistributed( x )
%ISDISTRIBUTED Summary of this function goes here
%   Detailed explanation goes here
isdistributedPaths = which('isdistributed','-ALL');
builtInIsdistributedPath = isdistributedPaths{end};
builtInIsdistributedPath(length(builtInIsdistributedPath)-15:length(builtInIsdistributedPath)) = [];
newIsdistributedPath = isdistributedPaths{1};
newIsdistributedPath(length(newIsdistributedPath)-15:length(newIsdistributedPath)) = [];
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    addpath(builtInIsdistributedPath);
    tf = isdistributed(x);
    addpath(newIsdistributedPath);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    tf = false;
end
