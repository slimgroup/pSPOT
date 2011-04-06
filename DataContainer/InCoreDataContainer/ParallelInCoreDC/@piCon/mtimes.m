function y = mtimes(A,B,swp)
% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = B;
    B = A;
    A = tmp;
    clear('tmp');
end

% Multiply
if ~isa(A,'dataContainer') % Right multiply
    y = piCon(A*double(B));
    
    % Extract collapsed dimensions
    y.imdims = { y.imdims{1} B.imdims{2} };
            
elseif ~isa(B,'dataContainer') % Left multiply
    y = piCon(double(A)*B);
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{1} y.imdims{2} };
    
else % Both data containers
    y = piCon(double(A)*double(B));
    
    % Extract collapsed dimensions
    y.imdims = { A.imdims{1} B.imdims{2} };
end

end % mtimes