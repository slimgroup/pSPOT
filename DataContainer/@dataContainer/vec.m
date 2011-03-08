function obj = vec(obj)

if ~obj.reallyveced
    
    tmp = obj;
    tmp.dims = [prod(tmp.dims) 1];
    tmp.veced = true;
    obj = tmp;
    clear('tmp');
end