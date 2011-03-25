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
            
            % Set default complexity
            if nargin == 1
                complex = false;
            end
            
            % Setup the file paths
            header = [f filesep 'FileProperties.xml'];
            trace  = [f filesep 'TraceFile'];
            dcpath = [f '.odcf' filesep];
            
            % Create datacon directory
            mkdir(dcpath);
            
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
            
            % Check for precision
            switch prec
                case 'float'
                    bytesize = 4;
                case 'double'
                    bytesize = 8;
                otherwise
                    error('Precision type not supported');
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
            dims = str2num([dims{:}]);
            
            % Close file
            fclose(fh);
            
            % Remove last dimension if it is single
            while(dims(end) == 1)
                dims(end) = [];
            end
            
            % Preallocate file
            if complex
                dims(1) = dims(1)/2;
                imagpath = [dcpath 'imag'];
                outCon.preallocateFile(imagpath, prod(dims));
            end
            realpath = [dcpath 'real'];
            outCon.preallocateFile(realpath, prod(dims));
            
            % Open trace file and Loop over last dimension
            ftrace    = fopen(trace);
            slicesize = prod(dims(1:end-1));
            for i = 1:dims(end)
                % Calculate offset and reset file
                offset  = slicesize*bytesize*(i-1);
                moffset = slicesize*8*(i-1);
                
                if complex
                    % Complex offset
                    fseek(ftrace,offset*2,-1);
                    
                    % Read real
                    x = fread(ftrace,slicesize,prec,bytesize);
                    
                    % Write real
                    M = memmapfile(realpath,'format',{'double',[slicesize 1],'A'},...
                        'offset',moffset, 'writable',true);
                    M.data(1).A = x;
                    
                    % Read imag
                    fseek(ftrace,offset*2 + bytesize,-1);
                    x = fread(ftrace,slicesize,prec,bytesize);
                    
                    % Write imag
                    M = memmapfile(imagpath,'format',{'double',[slicesize 1],'A'},...
                        'offset',moffset, 'writable',true);
                    M.data(1).A = x;                    
                    
                else
                    % Real offset
                    fseek(ftrace,offset,-1);
                    
                    % Read real
                    x = fread(ftrace,slicesize,prec);
                    
                    % Write real
                    M = memmapfile(realpath,'format',{'double',[slicesize 1],'A'},...
                        'offset',moffset, 'writable',true);
                    M.data(1).A = x;
                end
                delete(M);
            end
            fclose(ftrace);
            
            % Create header file
            fhead = fopen([dcpath 'head'],'w');
            
        end % javaSeisRead
        
        % allocFile - File preallocation function
        function preallocateFile(filename,numelements)
            
            % Check for numelements
            assert(isscalar(numelements) && isnumeric(numelements),...
                'Number of elements must be scalar and numeric.')
            
            % Allocate file
            allocFile(filename,numelements,8);
        end % allocFile
        
    end % Static methods
    
end % classdef
























