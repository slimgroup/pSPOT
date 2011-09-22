classdef (InferiorClasses = {?opSpot}) oppSpot < opSpot
   %oppSpot pSpot operator super class.
   %
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties
      gather = 0;
      weights;
      ddistscheme = []; % Distribution scheme of operator domain
      rdistscheme = []; % Distribution scheme of operator range
   end %properties
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Public methods
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      function op = oppSpot(type,m,n)
         %oppSpot  Constructor.
         op = op@opSpot(type,m,n);
         op.ddistscheme = n;
         op.rdistscheme = m;
      end
      
   end %methods - public
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Protected methods
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
%       function x = applyMultiply(op,x,mode)
%          op.counter.plus1(mode);
%          x = op.multiply(x,mode);
%          
%          
%          
%       end
      
      function x = applyDivide(op,x,mode)
         x = divide(op,x,mode);
      end
      
      % Signature of external protected functions
      x = divide(op,x,mode);
      
   end %methods - protected
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Abstract methods -- must be implemented by subclass.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Abstract, Access = protected )
       x = multiply(op,x,mode)
   end % methods - abstract
   
end % classdef
