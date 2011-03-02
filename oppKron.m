% classdef oppKron < oppKron2Lo
%     %OPPKRON    Kronecker tensor product to act on a distributed vector.
%     %   we need an oppKron that can be robust about the distribution
%     %   dimension.
%     %   oppKron(OP1,OP2,...,OPN,GATHER) Where OPi are Spot operators or
%     %   numeric matrices. Note that the last operator is applied to x then
%     %   the rest of the operators to the tranpose of the result, ans so x
%     %   should be of dimensions [cols(OP1),cols(OP2),...,cols(OPN)], and
%     %   vectorized after distribution so it is distributed along the last
%     %   dimension evenly.
%     %
%     %   Optional parameter gather specifies whether the output vector
%     %   should be gathered to the local lab.
%     %   GATHER = 0
%     
% end % classdef