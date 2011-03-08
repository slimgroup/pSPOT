function obj = vec(obj)

if ~obj.reallyveced
    
    tmp       = obj;
    d         = [tmp.dims{:}];
    tmp.dims  = {prod(d) 1};
    tmp.veced = true;
    obj       = tmp;
    clear('tmp');
end