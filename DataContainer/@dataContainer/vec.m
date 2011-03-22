function x = vec(x)
%VEC    Vectorization of the data container
%
%   vec(X) vectorizes dataContainer(X)
%
%   If X is distributed, and distributed along any dimensions other than
%   the last, it will be redistributed to the last dimension before
%   being vectorized. No communication happens in the reshaping stage.
%
%   See also: unvec, double, reshape

% Remove implicit
univec(x);

% Setup variables
if ~x.veced
    if x.isdist
        ddim = x.codist.Dimension;
        if ddim ~= length(size(x))
            warning('dataCon:RedistributingX',...
                'Redistributing x to the last dimension');
            x = distriCon(x);
        end
    end
    x = reshape(1,x,[numel(x.data) 1]);
end
