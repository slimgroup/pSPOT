classdef oppDictionary < oppSpot
    %OPPDICTIONARY  Dictionary of concatenated operators that acts in
    %               parallel
    %
    %   D = oppDictionary(OP1,OP2,...OPn,GATHER) creates a dictionary operator
    %   consisting of the concatenation of all operators. The optional last
    %   parameter, gather, specifies whether to gather the final result to a
    %   local variable instead of a distributed vector, by default this is set
    %   to 0.
    %
    %   Each of the operators is multiplied with the corresponding part of x on
    %   seperate labs.
    %
    %   *Note - only spot operators can be used in an oppDictionary, pSpot
    %   operators act in parallel already, and cannot be used with
    %   oppDictionary.
    %
    %   **Note - As of now the operators will be distributed according to 
    %   the Matlab default codistribution scheme. The x vector should be 
    %   distributed with this in mind.
    %   You could check the distribution scheme of your oppDictionary
    %   object by using this code:
    %
    %       child = A.children;
    %       spmd, child, end;
    %           
    %   for more info, type 'help codistributor1d'
    %
    %   See also oppBlockDiag, oppNumBlockDiag, oppStack
    
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
            assert( all(m == m(1)), 'Operator sizes are not consistent');
            n = sum(n);
            real = cellfun(@isreal,opList);
            cflag = ~all(real);
            linear = cellfun(@(p) logical(p.linear), opList);
            linear = all(linear);
            
            % Construct
            op = op@oppSpot('pDictionary', m(1), n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = distributed(opList);
            op.sweepflag= true;
            op.gather   = gather;
            op.precedence= 1;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            opchildren = gather(op.children);
            str = ['[',char(opchildren{1})];
            
            for ops=opchildren(2:end)
                str = [str, ', ', char(ops{1})];
            end
            
            str = [str, ']'];
            
            clear opchildren;
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = double(op)
            %OPPDICTIONARY.DOUBLE Distributed doubling of oppDictionary
            %    A = double(op) will apply double to each child operator
            %    in oppDictionary, and return a distributed dictionary
            %    of explicit operators.
            opchildren = op.children;
            childm = op.m;
            opn = op.n;
            spmd
                opchildren = getLocalPart(opchildren);
                childn = 0;
                partition = codistributed.zeros(1,numlabs);
                globalsize = [childm opn];
                A = zeros(childm,childn);
                if ~isempty(opchildren)
                    for i = 1:length(opchildren)
                        child = opchildren{i};
                        childn = childn + child.n;
                    end
                    A = zeros(childm,childn);
                    partition(labindex) = childn;
                    k = 0;
                    for i = 1:length(opchildren)
                        child = opchildren{i};
                        n = child.n;
                        A(:,k+1:k+n) = double(child);
                        k = k+n;
                    end
                end
                partition = gather(partition);
                codist = codistributor1d(2,partition,globalsize);
                A = codistributed.build(A,codist,'noCommunication');
            end % spmd
        end % double
        
    end % Methods
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            % Setting up class variables
            opchildren = op.children;
            opgather = op.gather;
            if ~isdistributed(x)
                error('X is not distributed');
            end
            
            spmd
                codist = codistributor1d(2,[],[1,length(opchildren)]);
                local_ops = codist.globalIndices(2);
                
                % Get sizes of the ops to go on the lab
                [M,N]  = cellfun(@size, opchildren);
                n = N(local_ops);
                sN = cumsum(N);
                
                % Multiplication
                opchildren = getLocalPart(opchildren);
                x = getLocalPart(x);
                if ~isempty(opchildren)
                    B = opDictionary(opchildren{:});
                    if mode == 1
                        if B.n ~= size(x,1)
                            error('Local x size mismatch. Probably wrong distribution');
                        end
                        y = B*x;
                    else
                        if B.m ~= size(x,1)
                            error('Local x size mismatch. Probably wrong distribution');
                        end
                        y = B'*x;
                    end
                end
                
                % Summing the results
                if mode == 1
                    if labindex == 1   %sum all results on lab 1
                        y = global_sum(y);
                    else
                        labSend(y,1);
                    end
                    
                    if ~opgather
                        y = codistributed(y,1,codistributor1d());
                    end
                else % mode 2
                    if opgather
                        y = gcat(y,1);   %concatenate results
                    else
                        part = codistributed.build( sum(n), ...
                            codistributor1d( 2, ones(1,numlabs), [1 numlabs]) );
                        codist = codistributor1d( 1, part, [sN(end) 1]);
                        y = codistributed.build( y, codist );
                    end
                end % summing
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