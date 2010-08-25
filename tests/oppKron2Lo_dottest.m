function status = oppKron2Lo_dottest(k,verbose)
%dottest  Apply the "dot-test" to an operator using distributed matrices.
%
%   dottest(A) generates random vectors X and Y (from the domain and range
%   of A), and verifies that (A*X)'*Y = X'*(A'*Y) within a tolerance of
%   1E-10. This can help detect errors in the operator; it canot be used to
%   guarantee correctness. The function returns false when the test
%   succeeded and true if it failed.
%
%   dottest(A,K) applies the dottest K times. (Default K is 100).
%
%   dottest(A,K,VERBOSE), with VERBOSE ~= 0, will print some
%   diagnostic output.

%   Copyright 2008-2009, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spot

% Set default parameters
if nargin < 1, k = 100; end
if nargin < 2, verbose = 0; end

% Initialize
tol    = sqrt(eps);
err    = 0;
kpass  = 0;

m = randi(100); o = randi(50); p = randi(50);
A = opDCT(m);
B = opCurvelet(o,p);
M = oppKron2Lo(A,B);

if verbose
   fprintf('Performing dot test on operator: %s\n', char(M));
end

for i=1:k
    x   = distributed.randn(B.n,A.n);
    y   = distributed.randn(B.m,A.m);
    x   = x(:);
    y   = y(:);
    z1  = (M*x)' * y;
    z2  = x' * (M'*y);
    err = max( err, abs(z1 - z2) );
    if err < tol,
       kpass = kpass + 1;
    end
end

if verbose
   if kpass < k
      fprintf('FAILED on %d out of %d tests\n', k-kpass, k);
      fprintf('%8s maximum absolute difference of %13.9e\n','',err);
   else
      fprintf('PASSED!\n');
   end
end

status = kpass == 0;
