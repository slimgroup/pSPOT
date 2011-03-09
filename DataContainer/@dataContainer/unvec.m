function y = unvec(x)
%UNVEC  Reshapes the data container into its pre-vectorized form
%
%   unvec(A) reshapes data container A using its original dimensions at the
%   time of construction. If A is distributed, this function will also
%   redistribute A back into the original distribution dimension. The
%   original distribution partition however, may not be conserved.
%
%   See also: vec, reshape, double

y = x;
if x.reallyveced    
    y             = reshape(x,x.odims{:});
    y.veced       = false;
    y.reallyveced = false;
    
else
    y.veced = false;
end