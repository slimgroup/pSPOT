function y = double(obj)

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
    
