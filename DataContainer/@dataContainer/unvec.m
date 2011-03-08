function obj = unvec(obj)

if obj.reallyveced
    
    tmp = obj;
    tmp.Data = reshape(tmp,tmp.odims);
    tmp.dims = tmp.odims;
    tmp.veced = true;
    obj = tmp;
    clear('tmp');    
end