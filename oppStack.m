function op = oppStack( varargin )
%OPPSTACK Summary of this function goes here
%   Detailed explanation goes here

if pSPOT.utils.hasDistComp % parallel installed
    op = oppStack_internal(varargin{:});
else
    warning(['Parallel Computing Toolbox not found,'...
        ' using serial version of operator']);
    % Remove gather
    if isscalar(varargin{end}) && ~isa(varargin{end}, 'opSpot')
        varargin(end) = [];
    end
    
    % Check for weights
    nargs = length(varargin);

    if isnumeric(varargin{1}) % weights
        weights = varargin{1};
        weights = weights(:);

        if nargs == 2 % Repeating ops                    
            if spot.utils.isposintscalar(varargin{1}) % repeat N times
                weights = ones(weights,1);                        
            end % Else: Repeating as many times as there are weights

            for i = 3:length(weights)+1
                varargin{i} = varargin{2};
            end
            varargin(1) = []; % delete weights  
        end
    end
    op = opStack(varargin{:});
end
