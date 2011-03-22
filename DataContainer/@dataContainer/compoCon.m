function x = compoCon(x)

assert(x.isdist,'x must be distributed')

data = x.data;
spmd
   data  = getLocalPart(data); 
end
x.data   = data;
x.compos = true;