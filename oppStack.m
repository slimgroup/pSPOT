classdef oppStack < oppSpot
%OPPSTACK  Stack of vertically concatenated operators.
%
%   oppStack(WEIGHTS, OP1, OP2, ...OPn,GATHER) creates a stacked operator
%   consisting of the vertical concatenation of all operators. When applied
%   the operators are divided amongst the labs and applied locally on each
%   lab. The optional last parameter, gather, specifies whether to gather
%   the final result to a local variable instead of a distributed vector,
%   by default this is set to 0. 
%
%               [WEIGHT1*OP1
%                WEIGHT2*OP2
%                   ...
%                WEIGHTn*OPn]
%
%   If the same weight is to be applied to each operator, set
%   WEIGHTS to a scalar. When WEIGHTS is empty [], it is set to
%   one. The WEIGHT parameter can be omitted as long as OP1 is not
%   a vector of length (n-1); in which case there is no way to
%   decide whether it is a weight vector or operator.
%
%   See also oppDictionary, opFoG, opSum.
%   
%   *Note - only spot operators can be used in an oppStack, pSpot
%   operators act in parallel already, and cannot be used with
%   oppDictionary.

%   Nameet Kumar - Oct 2010


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        weights;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppStack(varargin)
        
            % Checks weights parameter
            if ~isnumeric(varargin{1})
                weights = ones(nargin,1);
            else
                weights = varargin{1};
                if isempty(weights), weights = 1; end;
                [m,n] = size(weights);
                if (((m == 1) && (n == nargin-1)) || ...
                        ((n == 1) && (m == nargin-1)) || ...
                        ((m == 1) && (n == 1)))
                    weights = ones(nargin-1,1).*weights(:);
                    varargin = varargin(2:end);
                else
                    weights = ones(nargin,1);
                end
            end
                        
            % Check for gather parameter
            if isscalar( varargin{end} ) && any(varargin{end} == [0 1])
                gather = varargin{end};   
                varargin = varargin(1:end-1);
            else
                gather = 0;
            end
            
            % Check for empty operators and remove them
            ops = ~cellfun(@isempty,varargin);
            assert(any(ops),'At least one operator must be specified.');
            arrayfun(@(ind) warning('input "%d" is empty',ind), find(~ops));
            opList = varargin(ops);
            
            %Check for pSpot operators
            ops = cellfun(@(p) isa(p,'oppSpot'), opList);
            assert(~any(ops),' oppSpot operators are not supported');
            
            % Convert all Arguments to operators
            ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
            opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));
            
                        
            % Check consistency and complexity
            [m,n] = cellfun(@size,opList);
            assert( all(n == n(1)), 'Operator sizes are not consistant');
            m = sum(m);
            real = cellfun(@isreal,opList);
            cflag = ~all(real);
            linear = cellfun(@(p) logical(p.linear), opList);
            linear = all(linear);
            
            % Construct
            op = op@oppSpot('pStack', m, n(1));
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.sweepflag= true;
            op.gather   = gather;
            op.precedence= 1;
            op.weights  = weights;
            
        end %Constructor
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Display
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function str = char(op)
          % Initialize
          str = ['[',char(op.children{1})];
       
          for ops=op.children(2:end)
             str = [str, '; ', char(ops{1})];         
          end
          
          str = [str, ']'];
       end % Display
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Double
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function A = double(op)
          A = zeros(op.m,op.n);
          k = 0;
          for child=op.children
             m = child.m;
             A(k+1:k+m,:) = double(child);
             k = k + m;
          end
       end % double

    end % Methods
       
        
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          
          spmd
              % create a codist, and get global indices to see how many ops
              % on each lab
              codist = codistributor1d(2,[],[1,length(op.children)]);
              local_ops = codist.globalIndices(2);
              
              % Get sizes of the ops to go on the lab
              [m,n] = cellfun(@size, op.children(local_ops));
              M = gcat(m);
              sM = cumsum(M);

              if mode == 1
                  
                 y = zeros(sum(m),1);   %preallocate local results
                 for ops = local_ops
                    ind = [ -M(ops)+1, 0] + sM(ops)-sM(local_ops(1))+m(1);
                    y( ind(1):ind(2) ) = op.weights(ops) * applyMultiply(...
                        op.children{ops}, x, 1 );
                 end
                 
                 if op.gather
                     y = gcat(y,1);   %concatenate results
                 else
                     part = codistributed.build( sum(m), ...
                        codistributor1d( 2, ones(1,numlabs), [1 numlabs]) );
                     codist = codistributor1d( 1, part, [sM(end) 1]);
                     y = codistributed.build( y, codist );
                 end
                 
              else  % mode 2
                  
                 y = zeros(n(1),1);
                 for ops = local_ops
                    ind = [ -M(ops)+1, 0] + sM(ops);
                    xd = op.weights(ops) * x( ind(1):ind(2) );
                    y = y + applyMultiply( op.children{ops}, xd, 2 );
                 end
                 
                 if labindex == 1   %sum all results on lab 1
                     y = global_sum(y);
                 else
                     labSend(y,1);
                 end
                 
                 if ~op.gather
                     y = codistributed(y,1,codistributor1d());
                 end
                 
              end
              
          end %spmd
          if op.gather, y = y{1}; end    %if we gathered, the data is on lab1
              
       end % Multiply          

    end % Protected Methods
   
end % Classdef

function y = global_sum(y)
% loop through labs checking if its ready to send data, and recieve if it
% is

    labs = 2:numlabs;
    while ~isempty(labs)
        if( labProbe(labs(1)) )
            y = y + labReceive(labs(1));
            labs(1) = [];
        else
            labs = circshift( labs, [0 -1] );
        end
    end
    
end