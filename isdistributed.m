function tf = isdistributed( x )
%ISDISTRIBUTED Summary of this function goes here
%   Detailed explanation goes here
pSpotPath = pSPOT.path;
% Grab names of installed toolboxes
v = ver;
[inst{1:length(v)}] = deal(v.Name);

if ismember('Parallel Computing Toolbox', inst) % parallel installed
    rmpath(pSpotPath)
    %tf = builtin(isdistributed, x);
    tf = isdistributed(x);
    addpath(pSpotPath)
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of function']);
    tf = false;
end
