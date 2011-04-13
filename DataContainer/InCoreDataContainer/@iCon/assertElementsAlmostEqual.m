function assertElementsAlmostEqual(varargin)

varargin = cellfun(@(p) InCoreDataContainer.stripicon(p),...
    varargin,'UniformOutput',false');
assertElementsAlmostEqual(varargin{:});