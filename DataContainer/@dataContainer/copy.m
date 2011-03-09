function copy(y,obj)
%COPY   Copy data container properties
%
%   copy(A,B) copies all object properties from data container B to data
%   container A except the current dimensions and the explicit data.

assert(isa(y,'dataContainer') && isa(obj,'dataContainer'),...
    'Object must be a dataContainer')

y.odims       = obj.odims;
y.veced       = obj.veced;
y.reallyveced = obj.reallyveced;
y.isdist      = obj.isdist;
y.ddims       = obj.ddims;