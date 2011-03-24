classdef mdCon < dataContainer
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
        dirname  = '';
        format   = 'double';
        istemp   = false;
        odims    = [];
        
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
            x = x@dataContainer(zeros(dims{:}));
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
        
        
    end % Static methods
    
end % classdef