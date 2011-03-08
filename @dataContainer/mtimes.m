function y = mtimes(A,D)
    if ~isa(A,'dataContainer')
        y = A*D.Data;
        
    elseif ~isa(D,'dataContainer')
        y = A.Data*D;
        
    else % Both data containers
        y = A.Data*D.Data;
    end