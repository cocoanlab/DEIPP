function resize_axes(h_axes, sz_axes)

axes_pos = get(h_axes, 'Position');
screen_pos = get(0, 'ScreenSize');
max_iter = 10;

for iter_i = 1:max_iter
    fig_newpos = [1, 1, ...
        sz_axes(1) / axes_pos(3), ...
        sz_axes(2) / axes_pos(4)];
    if any(screen_pos(3:4) < fig_newpos(3:4))
        error('The required size of figure is too large!');
    end
    set(h_axes.Parent, 'Position', fig_newpos);
    axes_newpos = get(h_axes, 'Position');
    if axes_newpos == axes_pos
        break
    else
        axes_pos = axes_newpos;
    end
end

end