function y = blkdiag(varargin)
%BLKDIAG   Block-diagonal concatenation of operator input arguments.
%
%   B = blkdiag(OP1,OP2,...) produces the block-diagonal operator
%
%            [ OP1                ]
%        B = [      OP2           ]
%            [           ...      ]
%            [                OPN ]
%

varargin = cellfun(@double,varargin,'UniformOutput',false);
y = dcInCore(blkdiag(varargin{:}));
