function y = dzeros(A,varargin)
%DZEROS Distributed zero vector in operator domain
%
%   y = dzeros(A) generates a zero vector with the size of the operator
%   domain, distributed accordingly to the operators needs so that A*y is
%   a valid operation.
%
%   y = drandn(A,NCOLS) generates NCOLS vectors

if isempty(varargin)
    ncols = 1;
else
    ncols = varargin{1};
end

scheme = A.opsn;
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
        ypart   = zeros(scheme(labindex),ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme(labindex);
        ycodist = codistributor1d(1,ygpart,[sum(scheme),ncols]);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
else % Non-distributed
    y = zeros(A.n,ncols);
end
