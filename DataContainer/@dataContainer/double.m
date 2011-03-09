function y = double(x)
%DOUBLE     Returns the explicit representation of the data container
%
%   double(A) returns the underlying explicit representation of data
%   container A. If A is implicitly vectorized, the vectorization will
%   happen now, and stored explicitly in the data container. It uses
%   Matlab's default vectorization methods for this purpose.
%
%   Vectorization can be reversed by calling unvec on the data container
%   after this operation.
%
%   See also: unvec, vec, reshape

if x.reallyveced
    y = x.Data;
else
    if x.veced
        data          = x.Data;
        x.Data        = data(:);
        x.reallyveced = true;
        x.ddims       = 1;
        y             = x.Data;
        
    else
        y = x.Data;
    end
end
    
