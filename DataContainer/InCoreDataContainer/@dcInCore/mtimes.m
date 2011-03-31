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
    y = dcInCore(A*double(D));
    
    % Extract collapsed dimensions
    y.imdims = { y.imdims{1} D.imdims{2} };
            
elseif ~isa(D,'dataContainer') % Left multiply
    y = dcInCore(double(A)*D);
    
else % Both data containers
    y = dcInCore(double(A)*double(D));
end

end % mtimes