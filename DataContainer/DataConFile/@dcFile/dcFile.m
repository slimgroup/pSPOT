classdef dcFile < dataContainer
%DCFILE     Out-of-core subclass of Data Container
%
%   x = dcFile(FILENAME,DIMENSIONS,PERMISSION,OUTPUTNAME) creates an 
%   out-of-core data container x with data in file FILENAME having 
%   dimensions DIMENSIONS. OUTPUTNAME specifies the filename in which x 
%   will be saved if modified.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
       fname   = ''; % Filename of file
       outname = ''; % Filename of output file
       odims   = []; % Out of core dimensions
       fdims   = []; % Out of core dimensions for output
       iSlice  = 0;  % Index of out-of-core slice currently loaded into 
                     % memory, 0 for unloaded.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function x = dcFile(filename,dims,access,output)
            
            % Setup and extract variables
            if nargin < 4
                output = strcat('output_',filename);
            end
            if nargin < 3
                access = 'r+';
            end
            if nargin < 2
                error('Dimensions must be defined');
            end
            
            % Test if file is valid / create file
            [fid,msg] = fopen(filename,access);
            assert(fid > 0, msg);
            fclose(fid);
                        
            % Construct data container and update dimensions
            x         = x@dataContainer([]);
            x.odims   = dims;
            x.perm    = 1:length(dims);
            x.fname   = filename;
            x.outname = output;
            clearHistory(x);
            setHistory(x);
                        
        end % constructor
        
    end % public methods
    
end % classdef