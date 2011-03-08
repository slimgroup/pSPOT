function y = double(obj)
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

if obj.reallyveced
    y = obj.Data;
else
    if obj.veced
        data = obj.Data;
        obj.Data = data(:);
        obj.reallyveced = true;
        y = obj.Data;
        
    else
        y = obj.Data;
    end
end
    
