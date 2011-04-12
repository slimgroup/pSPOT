function y = invvec(x)

redims = [x.imdims{:}];
while(redims(end) == 1 && length(redims) > 2) % Strip singleton dimensions
    redims(end) = [];
end
y      = reshape(x,redims);