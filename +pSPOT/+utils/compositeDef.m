function y = compositeDef(x)
%COMPOSITEDEF Default compositioning of pspot operators
%   
%   y = compositeDef(x) returns a Composite array with the list of ops, x
%   stored in the default distribution configuration given the current
%   matlabpool conditions.

y       = Composite();
globind = pSPOT.utils.defGlobInd(length(x));

for i=1:length(y)
    y{i} = x(globind{i});
end