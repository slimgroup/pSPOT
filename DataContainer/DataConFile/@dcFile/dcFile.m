classdef dcFile < dataContainer
%DCFILE     Out-of-core subclass of Data Container
%
%   x = dcFile(FILENAME,DIMENSIONS,OUTPUTNAME) creates an out-of-core data
%   container x with data in file FILENAME having dimensions DIMENSIONS.
%   OUTPUTNAME specifies the filename in which x will be saved if modified.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
       fname   = ''; % Filename of file
       outname = ''; % Filename of output file
       odims   = []; % Out of core dimension
       iSlice  = 0;  % Index of out-of-core slice currently loaded into 
                     % memory, 0 for unloaded.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function x = dcFile(filename,dims,output)
            
            if nargin == 2
                output = strcat('output_',filename);
            elseif nargin < 2
                error('Dimensions must be defined');
            end
            
            % Test if file is valid
            [fid,msg] = fopen(filename,'r+');
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