function y = reform(x)

y = reshape(x.imcoddims,x,x.imdims{:});