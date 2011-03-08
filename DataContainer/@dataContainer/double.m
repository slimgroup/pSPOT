function y = double(obj)

if obj.reallyveced
    y = obj.Data;
else
    if obj.veced
        tmp = obj;
        tmp.Data = tmp.Data(:);
        tmp.reallyveced = true;
        y = tmp.Data;
        clear('obj');
    else
        y = obj.Data;
    end
end
    
