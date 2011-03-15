function x = unDistriCon(x)

if x.isdist
   x.data   = gather(x.data);
   x.isdist = false;
end