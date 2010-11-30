classdef oppBlockDiag < oppSpot
%OPPBLOCKDIAG   Operator-diagonal operator in parallel.
%
%   B = oppBlockDiag(OP1, OP2,...,OPN,OVERLAP,GATHER) creates a compound
%   block operator with the input operators OP1, OP2,... on the diagonal of
%   B, e.g., B = DIAG([OP1 OP2 ... OPN]). When multiplying the operators
%   are distributed among the labs and multiplied locally on each. When
%   OVERLAP is a positive integer the blocks will be offset OVERLAP rows
%   relative to the previous operator, when OVERLAP is negative the
%   operators are offset by the absolute value of OVERLAP in columns. Note
%   that choosing OVERLAP larger than the operator size may cause the
%   matrix to become block antidiagonal. GATHER specifies whether to gather
%   the results to a local array or leave them distributed, default is 0.
%
%   B = opBlockDiag(WEIGHT,OP1,...,OPN,OVERLAP,GATHER) additionally
%   weights each block by the elements of the vector WEIGHT. If
%   only a single operator is given it is replicated as many times
%   as there are weights.
%
%   B = opBlockDiag(N,OP,OVERLAP,GATHER) similar as above with WEIGHT
%   equal to ones(N,1). This will cause operator OP to be repeated
%   N times.
%
%   B = opBlockDiag(WEIGHT,A,OVERLAP,GATHER) where A is a 3D numerical
%   matrix. This will split the slices of A along the diagonal of B.
%
%   See also opFoG, opKron, opDictionary.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        weights;
        overlap;
        offset = 0;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppBlockDiag(varargin)
            
            % Extract gather and overlap
            if isscalar( varargin{end} ) && nargin > 2
                arg1 = varargin{end-1};
                arg2 = varargin{end};
                if nargin > 2 && isscalar( arg1 )
                    overlap   = arg1;
                    gather    = arg2;
                    varargin(end-1:end) = [];
                else
                    overlap   = arg2;
                    gather    = 0;
                    varargin(end) = [];
                end
                clear arg1 arg2
            else
                gather = 0;
                overlap = 0;
            end
            if ~spot.utils.isposintscalar(abs(overlap)+1)
                error('Overlap must be an integer scalar.');
            end
            
            % Extract weight
            nargs = length(varargin);
            if ~isnumeric(varargin{1})
                weights = ones(nargs,1);
            else
                weights = varargin{1};
                if isempty(weights), weights = 1; end;
                [m,n] = size(weights);
                
                % Repeating Op
                if spot.utils.isposintscalar(weights)&& nargs == 2
                    weights = ones(weights,1);
                    varargin(1) = [];
                    
                    % Normal
                elseif (((m == 1) && (n == nargs-1)) || ...
                        ((n == 1) && (m == nargs-1)) || ...
                        ((m == 1) && (n == 1)))
                    weights = ones(nargs-1,1).*weights(:);
                    varargin(1) = [];
                    
                    % 3D input matrix
                elseif nargs == 1 && ndims(varargin{1}) == 3
                    nSlices = size(varargin{1},3);
                    weights = ones(size(varargin{1},3),1);
                elseif nargs == 2 && ndims(varargin{2}) == 3
                    nSlices = size(varargin{2},3);
                    if ((m == 1) && (n == nSlices)) || ...
                            (m == nSlices) && (n == 1)
                        weights = ones(nSlices,1).*weights(:);
                    else
                        weights = ones(nSlices,1);
                    end
                    varargin(1) = [];
                    
                else
                    weights = ones(nargs,1);
                end
            end
            
            % Check for empty operators and remove them
            ops = ~cellfun(@isempty,varargin);
            assert(any(ops),'At least one operator must be specified.');
            arrayfun(@(ind) warning('input "%d" is empty',ind), find(~ops));
            opList = varargin(ops);
            
            % Check for pSpot operators
            ops = cellfun(@(p) isa(p,'oppSpot'), opList);
            assert(~any(ops),' oppSpot operators are not supported');
            
            % Convert all Arguments to operators -except 3D matrix case
            if ~exist('nSlices')
                ops = cellfun(@(p) ~isa(p,'opSpot'), opList);
                opList(ops) = cellfun(@(p) {opMatrix(p)}, opList(ops));
            end
            
            % Check complexity and setup single Op cases
            if exist('nSlices')           % 3D matrix
                
                opA = opList{1};
                [m,n,p] = size(opA);
                m = m*p;    n = n*p;
                cflag = ~isreal(opA) || ~all(isreal(weights));
                linear = 1;
                
                for i = 1:nSlices
                    opList{i} = opA(:,:,i);
                end
                clear opA
                
            elseif length(opList) == 1    % Repeat 1 operator
                
                opA = opList{1};
                [m,n] = size(opA);
                m = m * length(weights);
                n = n * length(weights);
                cflag  = ~isreal(opA) || ~all(isreal(weights));
                linear = opA.linear;
                
                [opList{1:length(weights)}] = deal(opA);
                clear opA
                
            else
                
                [m,n] = cellfun(@size,opList);
                m = sum(m);     n = sum(n);
                real = cellfun(@isreal,opList);
                cflag = ~all(real);
                linear = cellfun(@(p) logical(p.linear), opList);
                linear = all(linear);
                
            end
            
            % Handle overlaps
            if overlap == 0
                offset = 0;
            elseif overlap < 0
                %overlap in columns
                [m,n] = cellfun(@size,opList);
                m = sum(m);
                column = [0, cumsum( n + overlap )];
                offset = min( column(1:end-1) );
                n = max( column(2:end) - overlap ) - offset;
            else
                %overlap in rows
                [m,n] = cellfun(@size,opList);
                n = sum(n);
                row = [0, cumsum( m - overlap )];
                offset = min( row(1:end-1) );
                m = max( row(2:end) + overlap ) - offset;
            end
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', m, n);
            op.cflag    = cflag;
            op.linear   = linear;
            op.children = opList;
            op.weights  = weights;
            op.overlap  = overlap;
            op.offset   = -offset;
            op.sweepflag= true;
            op.gather   = gather;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            str = 'BlockDiag(';
            
            if ~isnumeric( op.children{1} )
                for child = op.children
                    str = [str,char(child{1}),', '];
                end
            else
                [m,n] = cellfun(@size, op.children);
                for i=1:length(op.children)
                    str = [str,'Matrix(',int2str(m(i)),',', ...
                        int2str(n(i)),'), '];
                end
            end
            str = [str(1:end-2), ')'];
        end % Display
        
    end % Methods
    
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 1
                overlap = op.overlap;
            else
                overlap = -op.overlap;
            end
                        
            spmd
                % create a codist, and get global indices to see how many
                % ops on each lab
                codist = codistributor1d(2,[],[1,length(op.children)]);
                loc_ops = codist.globalIndices(2);

                % Get sizes of the ops to go on the lab
                if mode == 1
                    [M,N] = cellfun(@size, op.children);
                else
                    [N,M] = cellfun(@size, op.children);
                end
                m   = M(loc_ops);
                n   = N(loc_ops);
                
                % iM,iN - first indices of operators
                if overlap == 0
                    iM = [0, cumsum( M )] + 1;
                    iN = [0, cumsum( N )] + 1;
                elseif overlap < 0
                    iN = [0, cumsum( N + overlap )] + op.offset + 1;
                    iM = [0, cumsum( M )] + 1;
                elseif overlap > 0
                    iM = [0, cumsum( M - overlap )] + op.offset + 1;
                    iN = [0, cumsum( N )] + 1;
                end
                fM = iM + [M,0] - 1;    % last index of operators
                fN = iN + [N,0] - 1;
                loc_start = min( iM(loc_ops));
                loc_rows = max( fM(loc_ops)) - loc_start + 1;

                % do the actual calculations
                res = zeros( loc_rows, size(x,2) );
                for ops = loc_ops
                    st = iM(ops) - loc_start + 1;
                    en = st + M(ops) - 1;
                    weight  = op.weights(ops);
                    child   = op.children{ops};
                    if mode == 2, weight = conj(weight); end
                    
                    if isnumeric(child)
                        if mode == 2, child = child'; end
                        res( st:en,: ) = res( st:en,: ) + weight * ...
                            child * x( iN(ops):fN(ops),: );
                    else
                        res( st:en,: ) = res( st:en,: ) + weight * ...
                            applyMultiply( child, x( iN(ops):fN(ops),: )...
                            , mode);
                    end
                end

                if labindex == 1    %gather results on lab 1

                    labs = 2:numlabs;
                    %build cell array to keep track of where everyone's
                    %results go
                    part = cumsum( codist.Partition );
                    %indices for lab1 res
                    ind = { [ loc_start, loc_start + loc_rows - 1] }; 
                    for l = labs  %add in indices for results of other labs
                        ind{l} = [ iM(part(l-1)+1), fM(part(l)) ];
                    end

                    y = zeros( max( fM ) , size(x,2));
                    y( ind{1}(1): ind{1}(2),: ) = res;
                    % loop through labs and add their results when
                    % they're ready to send
                    while ~isempty(labs)
                        if( labProbe(labs(1)) )
                            li = ind{labs(1)};
                            y( li(1):li(2),: ) = y( li(1):li(2),: ) + ...
                                labReceive(labs(1));
                            labs(1) = [];
                        else
                            labs = circshift( labs, [0 -1] );
                        end
                    end                            

                else
                    labSend(res,1);
                    y = 0;
                end

                if ~op.gather       %distributed the result from lab 1
                    y = codistributed(y,1,codistributor1d());
                end
                
            end % spmd
            if op.gather, y = y{1}; end     %gathered data is on lab 1
            
        end % Multiply
        
    end % Protected Methods
    
end % Classdef