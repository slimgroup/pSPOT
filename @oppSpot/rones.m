function y = rones(A,varargin)
%RONES Distributed ones vector in operator range
%
%   y = rones(A) generates a ones vector with the size of the operator
%   range, distributed accordingly to the operators needs so that A'*y is
%   a valid operation.
%
%   y = rones(A,NCOLS) generates NCOLS vectors

if isempty(varargin)
    ncols = 1;
else
    ncols = varargin{1};
end

scheme = A.rdistscheme;
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
        ypart   = ones(scheme(labindex),ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme(labindex);
        ycodist = codistributor1d(1,ygpart,[sum(scheme),ncols]);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
else % Non-distributed
    y = ones(A.m,ncols);
end