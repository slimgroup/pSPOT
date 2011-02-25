classdef oppQ < oppSpot
    %OPPQ   The final missing piece in pSpot
    %
    %   Q = oppQ(S,F) where S is a composite of vectors, and F a
    %   function handle that takes in local parts of S and gives a local 
    %   part of the final answer.
    %   F has to take in S and x, local parts of the respective vectors.
    %
    %   For now it is assumed that the size of the final answer will be
    %   the same as the size of x.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        fun;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppQ(S,F)
            
            [m,n] = size(S);
            op = op@oppSpot('Q', m, n);
            op.children = {S};
            op.fun = F;
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str='???';
            
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes is overloaded so as to call multiplication on a
        % distributed array. This multiplication will do the expected 2D
        % transform on 'x'.
        % For the moment mtimes is only implemented for right
        % multiplication
        function y=mtimes(op,x)
            if ~isa(op,'oppQ')
                error('Left multiplication not taken in account')
            elseif isnumeric(x) || isdistributed(x)
                assert( isvector(x) , 'Please use vectorized matrix')
                y=op.multiply(x, 1 ); 
            elseif isa(x,'opSpot')
                y = opFoG(op,x);
            else
                error(['unsupported data type: ' class(x)]);
            end
        end
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            S = op.children{1};
            F = op.fun;
            spmd
                % Setup local parts
                local_x = getLocalPart(x);
                local_S = getLocalPart(S);
                
                % Preallocate y and apply function
                y = zeros(size(local_x));
                y = F(local_S,local_x);
                
                y = codistributed.build(y,getCodistributor(x),'noCommunication');
                
            end % spmd
                
        end % multiply
        
    end % Protected methods
    
end % classdef























