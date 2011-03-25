classdef outCon < dataContainer
    %MDCON  Memmapfile Data Container class
    %
    %   mdCon(FILENAME,SIZE,PARAM1,VALUE1,...)
    %
    %   Parameters:
    %   format - The precision of the data file. default 'double'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = protected)
        dirname = '';
        format  = 'double';
        istemp  = false;
        odims   = [];
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = mdCon(f,varargin)
            
            % Check arguments size
            assert(rem(length(varargin), 2) == 0,...
                'Param/value pairs must come in pairs.')
            
            % Read the header files and extract details
            
            
            % Extract sizes
            assert(isnumeric(s),'Sizes must be numeric')
            dims = num2cell(s);
            dims{end} = 0;
            
            % Construct in-core data container
            x = mdCon@dataContainer(zeros(dims{:}));
            x.odims = s;
            
            % Extract filename
            assert(ischar(f),'Directory name must be a string')
            x.dirname = f;
                        
            % Parse param-value pairs
            for i = 1:2:length(varargin)
                
                assert(ischar(varargin{i}),...
                        'Parameter at input %d must be a string.', i);
                
                fieldname = lower(varargin{i});
                switch fieldname
                    case {'format','ifilename','ofilename'}
                        x.(fieldname) = varargin{i+1};
                    otherwise
                        error('Parameter "%s" is unrecognized.', ...
                            varargin{i});
                end
                
            end                       
            
        end % constructor
        
        function delete(x)
            if x.istemp
               rmdir(x.dirname); 
            end
        end % delete
        
    end % public methods
    
    methods(Static)
        
        % JaveSeiz reader, to be written properly in the near future
        function javaSeisRead(f,complex)
            
            if nargin == 1
                complex = false;
            end
            
            % Setup the file paths
            header = [f filesep 'FileProperties.xml'];
            trace  = [f filesep 'TraceFile'];
            
            % Open up the header file to extract stuffs
            fh = fopen(header);
            
            % Read in precision
            loop = true;
            while(loop)
                prec = fgetl(fh);
                ind = strfind(prec,'TraceFormat');
                if ind
                    ind  = ind + 28;
                    prec = prec(ind:end);
                    prec = lower(prec(1:strfind(prec,' ') - 1));
                    loop = false;
                end
            end
            
            % Read in dimensions
            loop = true;
            while(loop)
                l = fgetl(fh);
                if strfind(l,'AxisLengths');
                    loop = false;
                end
            end
            
            dims = {};
            i = 1;
            while(1)
                dims{i} = fgetl(fh);
                if strfind(dims{i},' </par>');
                    dims(i) = [];
                    break;
                end
                i = i + 1;
            end
            dims = [dims{:}];
            
            % Close file
            fclose(fh);
            
        end
        
        function preallocateFile(filename,length)
           allocFile(filename,length)
        end
        
    end % Static methods
    
end % classdef
























