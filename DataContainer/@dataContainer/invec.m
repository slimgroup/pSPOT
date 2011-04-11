function y = invec(x)

redims = [x.imdims{:}];
while(redims(end) == 1) % Strip singleton dimensions
    redims(end) = [];
end
y      = reshape(x,redims);