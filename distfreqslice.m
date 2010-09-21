function Y = distfreqslice(X, mode)
%takes in a t,x,y,.. distributed nd-array and returns x,y,..,f freq
%slice volume
assert( isa(X,'distributed'), 'input must be distributed');
if mode == 1
    Y = fft(X);
    Y = DistPermute(Y,[2:ndims(Y) 1]);
    
else
    sizes = size(X);
    perm = [ 2:ndims(X), 1];
    spmd
        codist = getCodistributor(X);
        part = codist.Partition;
        dim = codist.Dimension;
        Y = fft(X);
        Y = getLocalPart(Y);
        Y = permute(Y,perm);
        sizes = circshift(sizes,[0 -1]);
        Y = codistributed.build(Y,codistributor1d(dim-1,part,sizes));
        Y = redistribute(Y,codistributor1d(ndims(Y)));
    end
end
end