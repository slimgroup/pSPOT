function y = getSlice(x,ind,prec,skippy)
%GETSLICE   Brings a slice of the file into the core
%
%   y = getSlice(X,ISLICE,[PRECISION],[SKIP]) returns the ith slice of the 
%   out-of-core data stored in dcFile X. Note that ISLICE represents the 
%   index of the last dimension of dcFile X.
%
%   After being read in once, the data would be in-core and could be
%   accessed as a Matlab matrix.

% Process input
precision = 'double';
skip      = 0;
if nargin > 2
    precision = prec;
elseif nargin > 3
    skip      = skippy; 
end

% Open file
[fid,msg] = fopen(x.fname,'r');
assert(fid > 0, msg);

% Setup variables
dims      = x.odims;
chunksize = prod(dims(1:end-1));
offset    = chunksize*(ind-1)*8;

% Point to offset
fseek(fid,offset,-1);

% Read data
y         = fread(fid,chunksize,precision,skip);
fclose(fid);

% Reshape
shapesize = num2cell(dims(1:end-1));
y         = reshape(y,shapesize{:});

% Setup variables
x.data    = y;
x.dims    = size(y);
x.iSlice  = ind;
setHistory(x);