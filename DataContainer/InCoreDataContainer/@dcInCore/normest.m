function [e,cnt] = normest(S,tol)
%NORMEST Estimate the matrix 2-norm.
%
%   normest(S) is an estimate of the 2-norm of the matrix S.
%
%   normest(S,tol) uses relative error tol instead of 1e-6.
%
%   [nrm,cnt] = normest(..) also gives the number of iterations used.
%
%   This function is a minor adaptation of Matlab's built-in NORMEST.
%
%   See also NORM, COND, RCOND, CONDEST.

if nargin == 1, tol = 1e-6; end

[e,cnt] = normest(S.data,tol);