classdef oppWindow2D < oppSpot
    %OPPWINDOW2D 2D Windowing on distributed data
    %   op = oppWindow2D(GLOBALSIZE,LABSHAPE) returns an operator that 
    %   takes in distributed array of size GLOBALSIZE and distribute it 2D
    %   according to LABSHAPE (a vector specifying the number of labs in
    %   the 1st and 2nd dimension respectively)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties  Properties  Properties  Properties  Properties          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        global_size;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods  Public Methods  Public Methods  Public Methods      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor  Constructor  Constructor  Constructor  Constructor %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppWindow2D(gsize,labshape,varargin)
            
            % Checks and assertions
            assert(matlabpool('size') == 4,...
                '4 and only 4 labs are supported as of now.');
            assert(all(labshape == [2 2]), 'LABSHAPE must be [2 2] for now');
            assert(all(size(gsize) == [1 2]), 'GLOBALSIZE must be a 1 x 2 vector');
            assert(mod(gsize(1),2) == 0, 'gsize dimensions must be even');
            assert(mod(gsize(2),2) == 0, 'gsize dimensions must be even');
            
            % Build operator
            op             = op@oppSpot('Window2D',gsize(1)/2, gsize(1));
            op.global_size = gsize;
            op.sweepflag   = true;
            op.cflag       = 1;
            op.linear      = 1;
            
        end % constructor
    end % Public methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected Methods  Protected Methods  Protected Methods             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply  Multiply  Multiply  Multiply  Multiply  Multiply      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 1
                % Checks and assertions
                if isdistributed(x)
                    spmd, cod = getCodistributor(x); end
                    cod = cod{1};
                    assert(cod.Dimension == 2,...
                        'x must be distributed in the second dimension');
                    assert(all(cod.Partition == cod.Partition(1)),...
                        'Partition size must be the same across all labs');
                    assert(all(cod.Cached.GlobalSize == op.global_size),...
                        'Global size of data must be consistent with operator');
                else
                    error('x must be distributed')
                end
                
                % Prepare variables
                m = op.global_size(1);
                n = op.global_size(2);
                
                % SPMD
                spmd
                    loc_x = getLocalPart(x);
                    
                    % assign lab partners and swap data
                    switch labindex
                        case 1
                            rec_x = labSendReceive(2,2,loc_x(m/2+1:end,:));
                            loc_x = [loc_x(1:m/2,:) rec_x];
                        case 2
                            rec_x = labSendReceive(1,1,loc_x(1:m/2,:));
                            loc_x = [rec_x loc_x(m/2+1:end,:)];
                        case 3
                            rec_x = labSendReceive(4,4,loc_x(m/2+1:end,:));
                            loc_x = [loc_x(1:m/2,:) rec_x];
                        case 4
                            rec_x = labSendReceive(3,3,loc_x(1:m/2,:));
                            loc_x = [rec_x loc_x(m/2+1:end,:)];
                    end
                    
                    % Rebuild codistributed
                    codist = codistributor1d(2,(n/2)*ones(1,numlabs),[m/2 2*n]);
                    y = codistributed.build(loc_x,codist,'noCommunication');
                end % spmd
                
            else % mode 2
                
            end % if mode == 1
            
        end % Multiply
        
    end % Protected methods
    
end % classdef