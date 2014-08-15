function y = distributed( x )
%DISTRIBUTED Summary of this function goes here
%   Detailed explanation goes here

% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    y = builtin(distributed, x);
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    y = x;
end
