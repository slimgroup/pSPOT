classdef oplWindow1Dtpr < opSpot
%OPONES   Operator equivalent to ones function.
%
%   oplWindow1Dtpr(N,P,H)
    
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
       function op = oplWindow1Dtpr(varargin)
	  assert(nargin==3,'lWindow1Dtpr: wrong # of arguments')
	  n = varargin{1};
	  p = varargin{2};
	  h = varargin{3};
          m = n + (p - 1) * 2 * h;
	  op = op@opSpot('lWindow1Dtpr',m,n);
	  op.p = p;
	  op.h = h;
       end % function oplWindow1Dtpr
       
    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
	      [ A ] = pSPOT.pWindow.funWindow1DtprFor(op.n,op.p,op.h);
              y = A * x;
           else
	      [ B ] = pSPOT.pWindow.funWindow1DtprBck(op.n,op.p,op.h);
              y = B * x;
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
