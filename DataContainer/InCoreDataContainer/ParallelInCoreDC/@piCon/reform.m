function y = reform(x)

y = reshape(length([x.imdims{:}]),x,x.imdims{:});