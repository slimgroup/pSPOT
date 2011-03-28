function y = mtimes(A,D,swp)
% unswap
if nargin == 3 && strcmp(swp,'swap')
    tmp = D;
    D = A;
    A = tmp;
    clear('tmp');
end

% Multiply
if ~isa(A,'dataContainer') % Right multiply
    y = dataContainer(A*double(D));
    % Find out which dimensions collapsed
    n = size(y,2);
    
    % Go back in history and find out when the reshape happened
    % Setup variables
    dimsHistory = D.history.dims;
    
    % Travel back in history
    for i = (length(dimsHistory) - 1):-1:1
        if numel(dimsHistory{i}) > numel(dimsHistory{i + 1})
            unreshaped_dimension = dimsHistory{i}; % Time when dimension is
            % last changed
            break;
        elseif numel(dimsHistory{i}) == numel(dimsHistory{i + 1}) &&...
                all(dimsHistory{i} ~= dimsHistory{i + 1})
            unreshaped_dimension = dimsHistory{i}; % Time when dimension is
            % last changed
            break;
        end
    end
    
    % Check if the cumulative product of the dimension matches A.n
    s = 1;
    for j = 1:length(unreshaped_dimension)
        s = s*unreshaped_dimension(j);
        if s == size(A,2), break; end
    end
    
    % Concatenate into correct output size
    output_structure = [size(A,1) unreshaped_dimension(j+1:end)];
    
    % Append for single dimension case
    if length(output_structure) == 1
        output_structure = [output_structure 1];
    end
    
    % Save y's current dimensions and append history
    tempdims = y.dims;
    y.dims   = output_structure;
    y.perm   = 1:length(y.dims);
    setHistory(y);
    y.dims = tempdims;
    y.perm   = 1:length(y.dims);
    setHistory(y);
        
elseif ~isa(D,'dataContainer') % Left multiply
    y = dataContainer(double(A)*D);
    
else % Both data containers
    y = dataContainer(double(A)*double(D));
end

end % mtimes