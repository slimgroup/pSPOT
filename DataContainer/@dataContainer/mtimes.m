function y = mtimes(A,D,swp)
if nargin < 3
    swp = 'meh';
end
% unswap
if strcmp(swp,'swap')
    tmp = D;
    D = A;
    A = tmp;
    clear('tmp');
end

if ~isa(A,'dataContainer') % Right multiply
    y = A*double(D);

elseif ~isa(D,'dataContainer') % Left multiply
    y = double(A)*D;

else % Both data containers
    y = double(A)*double(D);
end