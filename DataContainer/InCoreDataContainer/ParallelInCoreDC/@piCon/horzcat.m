function y = horzcat(varargin)
%HORZCAT  Horizontal concatenation.
%
%   [A B] is the horizonal concatenation of data containers A and B.
%
%   See also piCon.vertcat

varargin = cellfun(@double,varargin,'UniformOutput',false');
y = piCon(horzcat(varargin{:}));
