function assertEqual(varargin)

varargin = cellfun(@(p) InCoreDataContainer.stripicon(p),...
    varargin,'UniformOutput',false');
assertEqual(varargin{:});