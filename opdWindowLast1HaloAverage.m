classdef opdWindowLast1HaloAverage < oppSpot
%opdWindowLast1HaloAverage tapered windowing for CARP-CG method
%
%   opdWindowLast1HaloAverage(N,L,P,H)
%
%   ARGUMENTS:
%      N = length of the input vector
%      L = length of the last dimension of input vector
%      P = number of processors
%      H = half of the overlap's size
%
%   ATTRIBUTES:
%      xshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the input vector in every window
%      yshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the output vector in every window
%
%   Notes:
%       1. This is a parallel/distributed operator
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        p = 0;
    l = 0;
    h = 0;
    yshape = 0;
    xshape = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opdWindowLast1HaloAverage(varargin)
      assert(nargin==4,'dWindowLast1HaloAverage: wrong # of arguments')
      n = varargin{1};
      l = varargin{2};
      p = varargin{3};
      h = varargin{4};
      assert(mod(n,l)==0,'Fatal error: L is not a valid last dimension')
      f = n/l;
          [ m xs ys ] = pSPOT.pWindow.funWindowLast1HaloShape( l, p, h );
      op = op@oppSpot('dWindowLast1HaloAverage',f*m,f*m);
      op.p = p;
      op.l = l;
      op.h = h;
      op.yshape = ys;
      op.xshape = xs;
       end % function opdWindowLast1HaloAverage
       
       % utest ( skip dottest )
       function output = utest(op,k,verbose)
           warning('No tests defined')
           output = 'PASSED!';
       end % utest

    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
           y = pSPOT.pWindow.funWindowLast1HaloAverage(x,op.l,op.p,op.h);
           else
            error('No transpose/inverse defined')
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
