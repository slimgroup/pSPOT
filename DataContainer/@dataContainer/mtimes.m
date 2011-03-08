function y = mtimes(A,D,swp)
if nargin < 3
    swp = 'meh';
end

if strcmp(swp,'swap')
    tmp = D;
    D = A;
    A = tmp;
    clear('tmp');
end

if ~isa(A,'dataContainer')
    y = dataContainer(A*double(D));

elseif ~isa(D,'dataContainer')
    y = dataContainer(double(A)*D);

else % Both data containers
    y = dataContainer(double(A)*double(D));
end