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
    if length(D.imdims) > length(y.imdims)
        
        collapsed_dims = 1;
        ind = 2;
        for i = 1:length(D.imdims)
            collapsed_dims = collapsed_dims * D.imdims(i);
            if collapsed_dims == size(A,2)
                ind = i + 1;
               break; 
            end
        end
        y.imdims = [y.imdims(1) D.imdims(ind:end)];
    end
            
elseif ~isa(D,'dataContainer') % Left multiply
    y = dcInCore(double(A)*D);
    
else % Both data containers
    y = dcInCore(double(A)*double(D));
end

end % mtimes