function writeSlice(x,y,precision,skip)

% Process input
if nargin < 4
    skip = 0;
end
if nargin < 3
    precision = 'double';
end

% Open file
[fid,msg] = fopen(x.outname,'a');
assert(fid > 0, msg);

% Setup variables
dims      = x.fdims;

% Write data
y         = fwrite(fid,y,precision,skip);
fclose(fid);
