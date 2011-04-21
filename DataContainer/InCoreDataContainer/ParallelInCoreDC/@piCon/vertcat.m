function y = vertcat(varargin)
%VERTCAT  Vertical concatenation.
%
%   [A; B] is the vertical concatenation of the data containers A and B.
%
%   See also piCon.horzcat

varargin = cellfun(@double,varargin,'UniformOutput',false');
y = piCon(vertcat(varargin{:}));
