classdef oppDictionary < oppSpot
%OPPDICTIONARY  Dictionary of concatenated operators that acts in
%parallel
%
%   D = oppDictionary(OP1,OP2,...OPn,GATHER) creates a dictionary operator
%   consisting of the concatenation of all operators. The optional last
%   parameter, gather, specifies whether to gather the final result to a
%   local variable instead of a distributed vector, by default this is set
%   to 0. 
%
%   ex.
%       D = oppDictionary(A,B,C);
%       x = D.drandn;
%       D*x;
%   or
%       x = rand( A.n*B.n*C.n, 1);
%
%   Each of the operators is multiplied with the corresponding part of x on
%   seperate labs.
%   
%   *Note - only spot operators can be used in an oppDictionary, pSpot
%   operators act in parallel already, and cannot be used with
%   oppDictionary.

%   Nameet Kumar - Oct 2010

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDictionary(varargin)
        
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
            assert( all(m == m(1)), 'Operator sizes are not consistant');
            n = sum(n);
            real = cellfun(@isreal,opList);
            cflag = ~all(real);
            linear = cellfun(@(p) logical(p.linear), opList);
            linear = all(linear);
            
            % Construct
            op = op@oppSpot('pDictionary', m(1), n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.sweepflag= true;
            op.gather   = gather;
            op.precedence= 1;
            
        end %Constructor
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Display
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function str = char(op)
          % Initialize
          str = ['[',char(op.children{1})];
       
          for ops=op.children(2:end)
             str = [str, ', ', char(ops{1})];         
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
             n = child.n;
             A(:,k+1:k+n) = double(child);
             k = k + n;
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
              [M,N]  = cellfun(@size, op.children);
              n = N(local_ops);
              sN = cumsum(N);

              if mode == 1
                  
                 y = zeros(M(1),1);
                 for ops = local_ops
                    ind = [ -N(ops)+1, 0] + sN(ops);
                    y = y + applyMultiply( op.children{ops},...
                         x( ind(1):ind(2) ), 1 );
                 end
                 
                 if labindex == 1   %sum all results on lab 1
                     y = global_sum(y);
                 else
                     labSend(y,1);
                 end
                 
                 if ~op.gather
                     y = codistributed(y,1,codistributor1d());
                 end
                 
              else  % mode 2
                  
                 y = zeros( sum(n), 1);     % preallocate local results
                 for ops = local_ops
                    ind = [ -N(ops)+1, 0] + sN(ops)-sN(local_ops(1))+n(1);
                    y( ind(1):ind(2) ) = applyMultiply( op.children{ops},...
                         x, 2 );
                 end
                 
                 if op.gather
                     y = gcat(y,1);   %concatenate results
                 else
                     part = codistributed.build( sum(n), ...
                        codistributor1d( 2, ones(1,numlabs), [1 numlabs]) );
                     codist = codistributor1d( 1, part, [sN(end) 1]);
                     y = codistributed.build( y, codist );
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