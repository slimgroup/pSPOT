classdef mdCon
    %MDCON  Memmapfile Data Container class
    %
    %   mdCon(FILENAME,PARAM1,VALUE1,...)
    %
    %   Parameters:
    %   format - The precision of the data file. default 'double'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = protected)
        filename  = '';
        format    = 'double';
        ifilename = '';
        ofilename = '';
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = mdCon(f,varargin)
            
            % Check arguments size
            assert(rem(length(varargin), 2) == 0,...
                'Param/value pairs must come in pairs.')
            
            % Extract filename
            assert(ischar(f),'Filename must be a string')
            x.filename = f;
            
            % Setup default output filename
            x.ofilename = strcat('output_',f);
            
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
        
    end % public methods
    
end % classdef