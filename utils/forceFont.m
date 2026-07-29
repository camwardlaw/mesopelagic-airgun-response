function forceFont(fig, fontName)
% Force a font onto every text/axes/legend object in a figure, regardless
% of how/when each object was created.
    if nargin < 2, fontName = 'Times New Roman'; end
    set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
end
