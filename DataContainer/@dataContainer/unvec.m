function x = unvec(x)
%UNVEC  Reshapes the data container into its pre-vectorized form
%
%   unvec(A) reshapes data container A using its original dimensions at the
%   time of vectorization. It goes back in history and fetches the correct
%   dimensions and codistributor to do this. Hence, the number of elements
%   and the distribution after vectorization must be conserved or unvec 
%   will fail horribly.
%
%   See also: vec, reshape, double

% Setup variables
dimsHistory = x.history.dims;
codHistory  = x.history.cod;
time = [];
% Travel back in history
for i = length(dimsHistory):-1:1
    if length(dimsHistory{i}) >= 2 && dimsHistory{i}(2) ~= 1
        time = i;
        break;
    end
end

% Error if pre-vectorized form not found
if isempty(time), error('A pre-vectorized form did not show up in history');end

% Check for number of elements
odims = dimsHistory{time};
ot    = time - 1;
pt    = time + 1;
vcod  = codHistory{pt}; % Distribution post-vec
ocod  = codHistory{ot}; % Codistributor before redistribution
assert(all(vcod.Partition == x.codist.Partition) &&...
prod(odims) == numel(x.data),'History has been altered beyond salvation')

% Reshape
x = reshape(length(odims),x,odims);

% Redistribute
x = distriCon(x,ocod);