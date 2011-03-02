function test_suite = test_oppDistFun
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
o = 5;
A = distributed.randn(m,n,o);
S = distributed.randn(m,o);
x = distributed.randn(n,o);
F = @funfun;
Q = oppDistFun(A,S,F);
x = x(:);
y = Q*x;

end %


%% Ingenious code that will forever solve the problem of last-dimension
% indexing, from
% http://www.mathworks.com/matlabcentral/answers/846-multidimensional-colon-operator
% idx(1:ndims(A) - 1) = {':'};
% A(idx{:},t) = x;
%