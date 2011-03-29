function y = mtimes(A,D,swp)
% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = D;
    D = A;
    A = tmp;
    clear('tmp');
end

% Multiply
if ~isa(A,'dataContainer') % Right multiply
    y = dataContainer(A*double(D));
            
elseif ~isa(D,'dataContainer') % Left multiply
    y = dataContainer(double(A)*D);
    
else % Both data containers
    y = dataContainer(double(A)*double(D));
end

end % mtimes