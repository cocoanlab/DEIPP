function vis_grp_r(r, grpidx, grpcols)

if ~isempty(grpidx)
    
    % ROI -> network
    [grp_sort, grp_sortidx] = sort(grpidx);
    imagesc(r(grp_sortidx, grp_sortidx));
    set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2, 'Clipping', 'off');

    grp_lines_offset = repmat(round(size(r,1) * 0.01), 2, numel(unique(grpidx)));
    grp_lines_pos = [0, find(diff(grp_sort)>0)' + 1; find(diff(grp_sort)>0)', numel(grp_sort)];
    grp_lines = line([-grp_lines_offset grp_lines_pos], [grp_lines_pos grp_lines_offset + numel(grp_sort)]);
    grp_lines_cols = repmat(grpcols(unique(grpidx), :), 2, 1);
    for line_i = 1:numel(grp_lines)
        set(grp_lines(line_i), 'LineWidth', 3, 'Color', grp_lines_cols(line_i, :));
    end
    
else
    
    % network
    imagesc(r);
    set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2, 'Clipping', 'off');
    for grp_i = 1:size(r,1)
        rectangle('Position', [-0.25 -0.25+grp_i 0.5 0.5], 'Curvature', [1 1], ...
            'FaceColor', grpcols(grp_i,:), 'EdgeColor', 'none');
        rectangle('Position', [-0.25+grp_i 0.75+size(r,1) 0.5 0.5], 'Curvature', [1 1], ...
            'FaceColor', grpcols(grp_i,:), 'EdgeColor', 'none');
    end
    
end

end