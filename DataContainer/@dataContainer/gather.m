function x = gather(x)

x.data   = gather(x.data);
x.isdist = false;
x.codist = [];
setHistory(x);