function cleanBox(ax)
% Keep the rectangular border look but drop the tick marks MATLAB mirrors
% onto the top/right edges whenever Box is 'on' -- turn box off, then
% redraw the 4-sided outline as a plain line (no ticks attached).
    xl = xlim(ax);
    yl = ylim(ax);
    box(ax, 'off');
    plot(ax, [xl(1) xl(2) xl(2) xl(1) xl(1)], [yl(1) yl(1) yl(2) yl(2) yl(1)], ...
        '-k', 'LineWidth', 0.5, 'HandleVisibility', 'off', 'Clipping', 'off');
    xlim(ax, xl);
    ylim(ax, yl);
end
