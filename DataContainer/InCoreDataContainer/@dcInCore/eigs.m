function varargout = eigs(varargin)
%EIGS   Find a few eigenvalues and eigenvectors of an operator using ARPACK.
%
%   eigs(A) returns six of the largest eigenvalues of an operator.
%
%   This routine is simply a wrapper to Matlab's own EIGS routine, and
%   most of the argument-list variations described in Matlab's EIGS
%   documentation are also allowed here.  The usage is identical to
%   Matlab's default version, except that the first argument must be a
%   Spot operator, and only the largest eigenvalues are considered.
%
%   Supported usage includes
%
%   D = EIGS(A)
%   [V,D] = EIGS(A)
%   [V,D,FLAG] = EIGS(A)
%   [V,D,FLAG] = EIGS(A,K)
%   [V,D,FLAG] = EIGS(A,K,SIGMA)
%   [V,D,FLAG] = EIGS(A,K,SIGMA,OPTS)
%
%   Again, see Matlab's built-in EIGS for details on these calls.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

   varargin{1} = varargin{1}.data;
   varargout{:} = eigs(varargin{:});