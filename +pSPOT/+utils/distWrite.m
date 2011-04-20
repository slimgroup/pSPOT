function distWrite(name,x,varargin)

% Setup variables
filename  = name;
precision = 'double';
repeat    = inf;
offset    = 0;

% Preprocess input arguments
error(nargchk(1, nargin, nargin, 'struct'));

if rem(length(varargin), 2) ~= 0
    error('Param/value pairs must come in pairs.');
end

assert(ischar(name), 'filename must be a string')
assert(isdistributed(x), 'Data must be distributed')

% Parse param-value pairs
for i = 1:2:length(varargin)
    
    assert(ischar(varargin{i}),...
        'Parameter at input %d must be a string.', i);
    
    fieldname = lower(varargin{i});
    switch fieldname
        case {'offset', 'precision', 'repeat'}
            eval([fieldname ' = varargin{i+1};']);
        otherwise
            error('Parameter "%s" is unrecognized.', ...
                varargin{i});
    end
end

% Set bytesize
switch precision
    case 'single'
        bytesize = 4;
        x = single(x);
    case 'double'
        bytesize = 8;
    otherwise
        error('Unsupported precision');
end

% Preallocate File
pSPOT.utils.allocFile(filename,prod(size(x)),bytesize);

spmd
    % Setup local chunk size
    local_data      = getLocalPart(x);
    local_size      = size(local_data);
    
    % Setup offsets
    global_elements = codistributed.zeros(1,numlabs);
    global_elements(labindex) = prod(local_size);
    global_elements = gather(global_elements);
    elements_offset = sum(global_elements(1:labindex-1));
    local_offset    = offset + elements_offset*bytesize;
    
    if labindex == 1   % Lab 1 always gets to write first
        % Setup memmapfile
        M = memmapfile(filename,'format',{precision,local_size,'x'},...
            'offset',local_offset, 'writable', true,...
            'repeat',repeat);
        
        % Write local data
        M.data(1).x = local_data;
    end
    labBarrier; % Synchronize
    
    if labindex == 1
        labs = 2:numlabs;
        while ~isempty(labs)
            if( labProbe(labs(1)) )
                labReceive(labs(1)); % Get the ready signal
                labSend('Go ahead',labs(1)); % Send the go ahead signal
                labReceive(labs(1)); % Wait for the done signal
                labs(1) = [];
            else
                labs = circshift( labs, [0 -1] );
            end
        end
    else % Other labs
        labSend('Im Ready!',1); % Send ready signal
        labReceive(1); % Wait for the go ahead signal
        % Setup memmapfile
        M = memmapfile(filename,'format',{precision,local_size,'x'},...
            'offset',local_offset, 'writable', true,...
            'repeat',repeat);
        
        % Write local data
        M.data(1).x = local_data;
        labSend('Im Done!',1);
    end
    
    % Setup memmapfile
    M = memmapfile(filename,'format',{precision,local_size,'x'},...
        'offset',local_offset, 'writable', true,...
        'repeat',repeat);
    
    % Write local data
    M.data(1).x = local_data;
    
end % spmd
