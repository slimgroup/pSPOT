function assertFalse(varargin)

varargin = cellfun(@(p) InCoreDataContainer.stripicon(p),...
    varargin,'UniformOutput',false');
assertFalse(varargin{:});