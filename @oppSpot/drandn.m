function y = drandn(A,varargin)
%DRANDN Distributed random vector in operator domain
%
%   y = drandn(A) generates a random vector y with the size of the operator
%   domain, distributed accordingly to the operators needs so that A*y is
%   a valid operation.
%
%   y = drandn(A,NCOLS) generates NCOLS vectors

if isempty(varargin)
    ncols = 1;
else
    ncols = varargin{1};
end

scheme = A.ddistscheme;
if length(scheme) > 1 % Distributed
    % scheme elements preprocessing
    nlabs = matlabpool('size');
    if length(scheme) > nlabs % If more than matlabpool, group
        newsindex = pSPOT.utils.defaultDistribution(length(scheme));
        j = 0; % Sum up the scheme elements
        for k=1:length(newsindex)
            newscheme(k) = sum(scheme(j+1:j+newsindex(k)));
            j = j+newsindex(k);
        end
        scheme = newscheme;
    else % nlabs >= length(scheme)
        % Append zeros at the back to account for empty labs
        scheme(end+1:end+(nlabs - length(scheme))) = 0;
    end
    spmd
        ypart = randn(scheme(labindex),ncols);
        ycodist = codistributor1d(1);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
else % Non-distributed
    y = randn(A.n,ncols);
end

