function part = partitionGenerator(dsize)
%PARTITIONGENERATOR Generates the default partition size for your
%                   distribtued array
%
%   PARTITION = partitionGenerator(DSIZE) generates a vector containing
%   the default partition based on the size of the desired distribution 
%   dimension.
%
%   Note that this function must be called from without a spmd block.
%   
%   See codistributor1d.defaultPartition for more details

nlabs = matlabpool('size');
part  = zeros(1,nlabs);

for i = 1:rem(dsize,nlabs)
    part(i) = ceil(dsize/nlabs);
end

for j = rem(dsize,nlabs) + 1 : nlabs
    part(j) = floor(dsize/nlabs);
end