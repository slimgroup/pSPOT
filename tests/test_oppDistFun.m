function test_suite = test_oppDistFun
% Experimental stage of oppDistFun
%test_oppCompositeFun  Unit tests for the oppCompositeFun operator
initTestSuite;
end

function test_oppDistFun_fun
%% Testing oppDistFun with a fun
% Solving the problem of indexing over last dimension

m  = 500;
n  = 300;
o  = 5;
A1 = distributed.randn(m,n,o);
A2 = distributed.randn(m,o);
x  = distributed.randn(n,o);
F  = @pSPOT.test.funfun;
Q  = oppDistFun(A1,A2,F);
x  = x(:);
y = Q*x;

end %

function test_oppDistFun_numBlockDiag
%% Testing a oppDistFun version of oppNumBlockDiag


end



%% Ingenious code that will forever solve the problem of last-dimension
% indexing, from
% http://www.mathworks.com/matlabcentral/answers/846-multidimensional-colon-operator
% idx(1:ndims(A) - 1) = {':'};
% A(idx{:},t) = x;
%