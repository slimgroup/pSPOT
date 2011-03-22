function x = spmdMultiply(A,x,mode)
%SPMDMULTIPLY   Multiplication across all dimensions higher than 2 in spmd
%   *This function must be called from within an spmd block
%
%   y = spmdMultiply(A,X) multiplies A with the the first and
%   second dimensions (slice) of the n-Dimensional codistributed array x
%   across all dimensions 3 to N recursively. OP can be a Spot operator or
%   a numerical matrix.

% Setup local parts
xcod  = getCodistributor(x);
gsize = size(x);

if xcod.Dimension == 1
    % Dimensional conflict
    % Redistribute to second dimension then redistribute back
    ncod  = codistributor1d(2);
    x     = redistribute(x,ncod);
    x     = getLocalPart(x);
    
    % Multiply
    if mode == 1
        x = spot.utils.nDimsMultiply(A,x);
    else
        x = spot.utils.nDimsMultiply(A',x);
    end
    
    % Rebuild and redistribute
    gsize(1) = size(x,1);
    ncod     = codistributor1d(2,ncod.Partition,gsize);
    x        = codistributed.build(x,ncod,'noCommunication');
    x        = redistribute(x,codistributor1d(1));
else
    % No dimensional conflict, proceed with multiplication
    x     = getLocalPart(x);
    ddims = xcod.Dimension;
    
    % Multiply
    if mode == 1
        x = spot.utils.nDimsMultiply(A,x);
    else
        x = spot.utils.nDimsMultiply(A',x);
    end
    
    % Rebuild
    gsize(1)     = size(x,1);
    ncod         = codistributor1d(ddims,xcod.Partition,gsize);
    labBarrier;
    x            = codistributed.build(x,ncod,'noCommunication');
end