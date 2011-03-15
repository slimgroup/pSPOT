function x = spmdMultiply(A,x,mode)
%SPMDMULTIPLY   Multiplication across all dimensions higher than 2 in spmd
%   *This function must be called from within an spmd block
%
%   y = spmdMultiply(A,X) multiplies A with the the first and
%   second dimensions (slice) of the n-Dimensional codistributed array x 
%   across all dimensions 3 to N recursively. OP can be a Spot operator or 
%   a numerical matrix.

% Setup local parts
xcod = getCodistributor(x);
if xcod.Dimension == 1
    % Dimensional conflict
    % Redistribute to second dimension then redistribute back
    x     = redistribute(x,codistributor1d(2));
    x     = getLocalPart(x);
    part  = codistributed.zeros(1,numlabs);
    
    % Multiply
    if mode == 1
        x = spot.utils.nDimsMultiply(A,x);
    else
        x = spot.utils.nDimsMultiply(A',x);
    end
    
    % Rebuild and redistribute
    part(labindex) = size(x,2);
    gsize = size(x);
    gsize(2) = sum(gather(part));
    ncod  = codistributor1d(2,part,gsize);
    x     = codistributed.build(x,ncod,'noCommunication');
    x     = redistribute(x,codistributor1d(1));
else
    % No dimensional conflict, proceed with multiplication
    part  = codistributed.zeros(1,numlabs);
    x     = getLocalPart(x);
    
    % Multiply
    if mode == 1
        x = spot.utils.nDimsMultiply(A,x);
    else
        x = spot.utils.nDimsMultiply(A',x);
    end
    
    % Rebuild
    part(labindex) = size(x,xcod.Dimension);
    gsize = size(x);
    gsize(xcod.Dimension) = sum(gather(part));
    ncod  = codistributor1d(xcod.Dimension,part,gsize);
    x     = codistributed.build(x,ncod,'noCommunication');
end