function x = distRead(name,dimensions,varargin)

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
assert(isnumeric(dimensions), 'dimensions must be numeric')

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
    case 'double'
        bytesize = 8;
    otherwise
        error('Unsupported precision');
end    

spmd
    % Setup local chunk size
    global_codist   = codistributor1d(length(dimensions),[],dimensions);
    partition       = global_codist.Partition;
    local_size      = [dimensions(1:end-1) partition(labindex)];
    
    % Setup offsets
    global_part     = codistributed(1:dimensions(end));
    global_indices  = globalIndices(global_part,2,labindex);
    elements_offset = prod([dimensions(1:end-1) global_indices(1) - 1]);
    local_offset    = offset + elements_offset*bytesize;
    
    % Setup memmapfile
    M = memmapfile(filename,'format',{precision,local_size,'x'},...
        'offset',local_offset,'repeat',repeat);
    
    % Read local data
    local_data      = double(M.data(1).x);
    x = codistributed.build(local_data,global_codist,'noCommunication');
    
end % spmd