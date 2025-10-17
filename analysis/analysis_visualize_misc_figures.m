function analysis_visualize_misc_figures(varargin)

preptype = 'cocoan-preproc';
antype = 'cocoan-analysis';
bidsdir = fileparts(fileparts(mfilename('fullpath'))); % mfilename: bidsdir/code/~.m
if isempty(bidsdir); bidsdir = fileparts(pwd); end
parctype = 'indparc';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'preptype'}
                preptype = varargin{i+1};
            case {'antype'}
                antype = varargin{i+1};
            case {'bidsdir'}
                bidsdir = varargin{i+1};
            case {'parctype'}
                parctype = varargin{i+1};
        end
    end
end

addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

sj_num = 1;

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
ratdir = fullfile(andir, 'rating', sj_id);
preddir = fullfile(andir, 'predict_rating', sj_id);
miscfigdir = fullfile(andir, 'misc_figures');

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat');

testname = fullfile(preddir, sprintf('%s_predict_rating_test_%s.mat', sj_id, parctype));
if ~exist(testname, 'file'); return; end
load(testname, 'Yfit');

figbasename = fullfile(miscfigdir, 'misc_figures');

%%

cols = call_colors;
colnames = {'wb_cw_w', 'wb_ry', 'wb_blb', 'wb_rybcbl', 'wb_rybcbl_blbc', 'wb_rybcbl_rybl'};
for i = 1:numel(colnames)
    vis_colorbar(cols.fx_interp(cols.(colnames{i}), 1000));
    pagesetup(gcf);
    exportgraphics(gcf, sprintf('%s_colorbar_%s.pdf', figbasename, colnames{i}), 'BackgroundColor', 'none');
    close all;
end

%%

rng(2019);
dat(:,1) = zscore(lowpass(randn(500,1),0.01));
dat(:,2) = zscore(lowpass(randn(500,1),0.01));
dat(:,3) = dat(:,1) .* dat(:,2);
cols = [58 170 224
    249 126 103
    249 155 0] ./ 255;

for i = 1:3
    plot(dat(1:375,i), 'LineWidth', 4, 'Color', cols(i,:));
    line([0;375], [0;0], 'LineWidth', 4, 'Color', cols(i,:));
    axis off; box off;
    set(gcf, 'Position', [680   911   432    85]);
    exportgraphics(gcf, sprintf('%s_edge_timeseries_ex_%d.pdf', figbasename, i), 'BackgroundColor', 'none');
    close all;
end

%%

rating = rating_dat{2,5}(1:1102);

plot(rating, 'LineWidth', 4);
box off; axis off;
set(gcf, 'color', 'w', 'position', [680   756   918   164]);
pagesetup(gcf);
exportgraphics(gcf, sprintf('%s_rating_ex.pdf', figbasename), 'BackgroundColor', 'none');
close all;

%%

cols = [255	131 91
    232 172 70
    15 227 128
    185 127 254] ./ 255;
nbin = 10;

hold on;
ylim([min(rating)-0.05 max(rating)+0.05]);
rating_interp = interp1([1:numel(rating)].', rating, linspace(1,numel(rating),numel(rating)*1000).');
rating_dscrank = rating_interp;
[~, wh_sortbin] = sort(rating_dscrank, 'descend');
rating_dscrank(wh_sortbin) = 1:numel(rating_dscrank);
rating_quantile = sum(rating_dscrank > [-Inf quantile(rating_dscrank, (1:nbin-1)./nbin)], 2);
for bin_i = [1 2 3 10]
    rating_quantile_bin = double(rating_quantile == bin_i);
    rating_quantile_bin_ends = reshape(find(ismember(diff([0; rating_quantile_bin; 0]), [-1 1])), 2, []).';
    rating_quantile_bin_ends(diff(rating_quantile_bin_ends, [], 2) < 1000*1, :) = [];
    patch([rating_quantile_bin_ends(:,1) rating_quantile_bin_ends(:,2) rating_quantile_bin_ends(:,2) rating_quantile_bin_ends(:,1)].', ...
        repmat([min(rating)-0.05 min(rating)-0.05 min(rating)-0.03 min(rating)-0.03].', 1, size(rating_quantile_bin_ends,1)), cols(min(bin_i,4),:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    rating_quantile_bin(~rating_quantile_bin) = NaN;
    plot(1:numel(rating_interp), rating_interp .* double(rating_quantile_bin), 'Color', cols(min(bin_i,4),:), 'LineWidth', 4);
end
rating_quantile_bin = double(~ismember(rating_quantile, [1 2 3 10]));
rating_quantile_bin(~rating_quantile_bin) = NaN;
plot(1:numel(rating_interp), rating_interp .* double(rating_quantile_bin), 'Color', [0.8 0.8 0.8], 'LineWidth', 4);

box off; axis off;
set(gcf, 'color', 'w', 'position', [680   756   918   164]);
pagesetup(gcf);
exportgraphics(gcf, sprintf('%s_rating_ex_bin.pdf', figbasename), 'BackgroundColor', 'none');
close all;

%%

h = histfit(rating_interp);
[x, y] = deal(h(2).XData, h(2).YData);
[x, y] = deal(x(x <= max(rating_interp) & x > min(rating_interp)), ...
    y(x <= max(rating_interp) & x > min(rating_interp)));
xlim([min(x) max(x)]);
delete(h);
hold on;
for bin_i = [1 2 3 10]
    idx = find(cumsum(y) <= max(cumsum(y)) .* bin_i/nbin & cumsum(y) >= max(cumsum(y)) .* (bin_i-1)/nbin);
    area(x(idx), y(idx), 'FaceColor', cols(min(bin_i,4),:), 'EdgeColor', cols(min(bin_i,4),:));
end
for bin_i = 4:9
    idx = find(cumsum(y) <= max(cumsum(y)) .* bin_i/nbin & cumsum(y) >= max(cumsum(y)) .* (bin_i-1)/nbin);
    area(x(idx), y(idx), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
end
h_area = findobj(gca, 'Type', 'Area');
for i = 1:numel(h_area); h_area(i).BaseLine.Visible = 'off'; end
box off; axis off;
set(gcf, 'color', 'w', 'position', [680   756   738   162]);
pagesetup(gcf);
exportgraphics(gcf, sprintf('%s_rating_ex_binhist.pdf', figbasename), 'BackgroundColor', 'none');
close all;

%%

pred = Yfit{2,5}(1:1102);
pred(773:795) = pred(773:795) .* 0.6;
pred = smooth(pred);

plot(pred, 'LineWidth', 4, 'Color', [0.8500    0.3250    0.0980]);
box off; axis off;
set(gcf, 'color', 'w', 'position', [680   756   918   164]);
pagesetup(gcf);
exportgraphics(gcf, sprintf('%s_pred_ex.pdf', figbasename), 'BackgroundColor', 'none');
close all;

%%

rng(2020);
a = randn(50, 1);
b = a + randn(50, 1).*1.5;
a = (a - min(a)) ./ (max(a) - min(a));
b = (b - min(b)) ./ (max(b) - min(b));
scatter(a, b, 200, [0.7 0.7 0.7], 'filled');
line(xlim, glmfit(a,b).'*[ones(1,2); xlim], 'LineWidth', 6, 'color', [0.7 0.7 0.7]);
box off; axis off; resize_axes(gca, [300 300]);
exportgraphics(gcf, sprintf('%s_eval_ex_corr.pdf', figbasename), 'BackgroundColor', 'none');
close all;

h = roc_plot(b, a > quantile(a,0.5));
box off; axis off; resize_axes(gca, [300 300]);
set(h.line_handle(1), 'MarkerSize', 10, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7]);
set(h.line_handle(2), 'LineWidth', 6, 'Color', [0.7 0.7 0.7]);
exportgraphics(gcf, sprintf('%s_eval_ex_roc.pdf', figbasename), 'BackgroundColor', 'none');
close all;

end