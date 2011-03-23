%% Experimenting with memmapfiles

% Create/overwrite data file
x = randn(10,10);
f = fopen('A','w');
fwrite(f,x,'double');
fclose(f);

M = memmapfile('A','format',{'double' [10 10] 'x'},...
    'writable',true)
%%
spmd
   if labindex == 1
       M.data.x(:,1:5) = opDFT(10)*M.data.x(:,1:5)
   else
       M.data.x(:,6:10) = opDFT(10)*M.data.x(:,6:10)
   end
end