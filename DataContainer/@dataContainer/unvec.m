function obj = unvec(obj)
%UNVEC  Reshapes the data container into its pre-vectorized form
%
%   unvec(A) reshapes data container A using its original dimensions at the
%   time of construction. If A is distributed, this function will also
%   redistribute A back into the original distribution dimension. The
%   original distribution partition however, may not be conserved.
%
%   See also: vec, reshape, double

if obj.reallyveced    
    obj             = reshape(obj,obj.odims{:});
    obj.veced       = false;
    obj.reallyveced = false;
    
else
    obj.veced = false;
end