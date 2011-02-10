function [opList,m,n,cflag,linear] = stdchk(varargin)
%STDCHK    Performs the standard checking and data extraction routine
%          in pSpot.
%
%   Returns:
%   opList: List of Spot operators in a cell array
%   [m,n]: The array of the sizes of the operators
%   cflag: Complexity of the operators
%   linear: Linearity of the operators

% Check for empty operators and remove them
ops = ~cellfun(@isempty,varargin);
assert(any(ops),'At least one operator must be specified.');
arrayfun(@(ind) warning('input "%d" is empty',ind), find(~ops));
opList = varargin(ops);

% Check for pSpot operators
ops = cellfun(@(p) isa(p,'oppSpot'), opList);
assert(~any(ops),' oppSpot operators are not supported');

% Convert all non-Spot arguments to operators
ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));

% Check complexity and setup sizes
[m,n] = cellfun(@size,opList);
real = cellfun(@isreal,opList);
cflag = ~all(real);
linear = cellfun(@(p) logical(p.linear), opList);
linear = all(linear);

end