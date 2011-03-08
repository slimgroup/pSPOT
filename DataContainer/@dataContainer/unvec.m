function obj = unvec(obj)

if obj.reallyveced    
    obj             = reshape(obj,obj.odims{:});
    obj.veced       = false;
    obj.reallyveced = false;
    
else
    obj.veced = false;
end