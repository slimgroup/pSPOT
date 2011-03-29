function x = empty(varargin)

if numel(varargin{1}) > 1
    varargin = num2cell(varargin{1});
end
varargin( end + 1 ) = {0};
x = dcInCore(zeros(varargin{:}));