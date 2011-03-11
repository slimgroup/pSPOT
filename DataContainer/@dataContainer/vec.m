function x = vec(x)
%VEC    Vectorization of the data container
%
%   vec(A) vectorizes dataContainer(A)
%
%   See also: unvec, double, reshape

% Setup variables
if x.isdist
    ddim = x.codist.Dimension;
    if ddim ~= length(size(x))
        warning('Redistributing x to the last dimension');
        x = distriCon(x);
    end    
end
x = reshape(1,x,[numel(x.data) 1]);
