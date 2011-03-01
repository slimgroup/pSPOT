function test_suite = exp_oppDistFun
% Experimental stage of oppDistFun
%test_oppCompositeFun  Unit tests for the oppCompositeFun operator
initTestSuite;
end

function test_oppDistFun_fun
%% Testing oppDistFun with a fun
% Solving the problem of indexing over last dimension
nlabs = matlabpool('size');
m = 500;
n = 300;
A = distributed.randn(m,n,nlabs);
S = distributed.randn(m,nlabs);
x = distributed.randn(n,nlabs);
F = @funfun;
Q = oppDistFun(A,S,F);

end %

function varargout = funfun(varargin)
%% Function for testing
if nargin == 1 && varargin{1} == 0
    y = [500 300 1 1];
    return;
else % Multiply
    A = varargin{1};
    S = varargin{2};
    x = varargin{3};
    mode = varargin{4};
    
    if mode == 1
        y = A*x - S;
    else
        y = A'*x - conj(S);
    end
end % multiply

end % funfun

%% Ingenious code that will forever solve the problem of last-dimension
% indexing, from
% http://www.mathworks.com/matlabcentral/answers/846-multidimensional-colon-operator
% idx(1:ndims(A) - 1) = {':'};
% A(idx{:},t) = x;
%