function tf = isdistributed( x )
%ISDISTRIBUTED Summary of this function goes here
%   Detailed explanation goes here
slimapps = getenv('SLIM_APPS');
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    rmpath([slimapps '/tools/utilities/pSPOT'])
    %tf = builtin(isdistributed, x);
    tf = isdistributed(x);
    addpath([slimapps '/tools/utilities/pSPOT'])
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    tf = false;
end
