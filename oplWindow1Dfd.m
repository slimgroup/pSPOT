classdef oplWindow1Dfd < opSpot
%OPONES   Operator equivalent to ones function.
%
%   oplWindow1Dfd(N,P,H)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        p = 0;
	h = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oplWindow1Dfd(varargin)
	  assert(nargin==3,'lWindow1Dfd: wrong # of arguments')
	  n = varargin{1};
	  p = varargin{2};
	  h = varargin{3};
          m = n + (p - 1) * 2 * h;
	  op = op@opSpot('lWindow1Dfd',m,n);
	  op.p = p;
	  op.h = h;
       end % function oplWindow1Dfd
       
    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
	      [ A ] = pSPOT.pWindow.funWindow1DfdFor(op.n,op.p,op.h);
              y = A * x;
           else
	      [ B ] = pSPOT.pWindow.funWindow1DfdBck(op.n,op.p,op.h);
              y = B * x;
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
