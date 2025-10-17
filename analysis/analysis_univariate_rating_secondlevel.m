function analysis_univariate_rating_secondlevel(sj_num, spmdir, mricrogldir, varargin)

preptype = 'cocoan-preproc';
antype = 'cocoan-analysis';
bidsdir = fileparts(fileparts(mfilename('fullpath'))); % mfilename: bidsdir/code/~.m
if isempty(bidsdir); bidsdir = fileparts(pwd); end
parctype = 'indparc';
mricroglclutdir = strrep(mricrogldir, fullfile('MacOS', 'MRIcroGL'), fullfile('Resources', 'lut'));

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
            case {'mricroglclutdir'}
                mricroglclutdir = varargin{i+1};
        end
    end
end

addpath(genpath(spmdir));
addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
parcdir = fullfile(andir, 'parcellation', sj_id);
scandir = fullfile(andir, 'SCAN', sj_id);
univdir = fullfile(andir, 'univariate_rating', sj_id);

cols = call_colors;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

vertnetfile = fullfile(scandir, sprintf('%s_rawassn_minsize400_regularized_recolored_addSCAN.dscalar.nii', sj_id));
if ~exist(vertnetfile, 'file'); return; end
vertnetcf = ft_read_cifti_mod(vertnetfile);
vertnetorig = parc.id;
vertnetorig(~isnan(parc.id)) = parc.netorig(parc.id(~isnan(parc.id)));
vertnetorig(1:numel(vertnetcf.data)) = vertnetcf.data;
[~, vertnet] = ismember(vertnetorig, [0:1 1.5 2:22]);

firstlvlfile = fullfile(univdir, sprintf('%s_univariate_rating_firstlevel_b.dscalar.nii', sj_id));
if ~exist(firstlvlfile, 'file'); return; end
firstlvlcf = ft_read_cifti_mod(firstlvlfile);

figbasename = fullfile(univdir, sprintf('%s_univariate_rating_secondlevel', sj_id));

%% Second-level t-test

[~, t_p, ~, t_stats] = ttest(firstlvlcf.data.');
t_p = t_p.';
t_t = t_stats.tstat.';

[t_fdr_h, t_fdr_p_thr] = fdr_bh(t_p, 0.05);
t_fdr_h_pos = t_fdr_h & t_t > 0;
t_fdr_h_neg = t_fdr_h & t_t < 0;

fprintf('Total %d vertices (pos %d, neg %d) survived at FDR p < 0.05 (uncorrected p < %.2f).\n', ...
    sum([t_fdr_h t_fdr_h_pos t_fdr_h_neg]), round(t_fdr_p_thr,2));

switch sj_num
    case 1
        t_fdr_thr = [-8.5 max(t_t(t_fdr_h_neg)) min(t_t(t_fdr_h_pos)) 8.5];
    case 2
        t_fdr_thr = [-7.6 max(t_t(t_fdr_h_neg)) min(t_t(t_fdr_h_pos)) 7.6];
end

%% Visualize cortex

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.mapname = {'unthresholded t-value'};
imgname = sprintf('%s_t_unthr', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 %f -neg-user 0 %f -palette-name RY-BC-BL', ...
    imgname, imgname, t_fdr_thr([4 1])));

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.data(~t_fdr_h) = NaN;
tempcf.mapname = {'thresholded t-value'};
imgname = sprintf('%s_t_fdr05', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 %f -neg-user 0 %f -palette-name RY-BC-BL', ...
    imgname, imgname, t_fdr_thr([4 1])));

tempcf = parc.cf;
tempcf.data = vertnetorig;
tempcf.mapname = {'unthresholded vertex network label'};
imgname = sprintf('%s_net_unthr', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));

tempcf = parc.cf;
tempcf.data = vertnetorig;
tempcf.data(~t_fdr_h) = NaN;
tempcf.mapname = {'thresholded vertex network label'};
imgname = sprintf('%s_net_fdr05', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));

tempcf = parc.cf;
tempcf.data = parc.id .* 0 + 1;
tempcf.mapname = {'unthresholded vertex binary'};
imgname = sprintf('%s_binary_unthr', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 2 -palette-name Gray_Interp_Positive', imgname, imgname));

tempcf = parc.cf;
tempcf.data = parc.id .* 0 + 1;
tempcf.data(~t_fdr_h) = NaN;
tempcf.mapname = {'thresholded vertex binary'};
imgname = sprintf('%s_binary_fdr05', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 2 -palette-name Gray_Interp_Positive', imgname, imgname));

%% Visualize cerebellum

tmpdir = tempname(tempdir); mkdir(tmpdir);

imgname = fullfile(tmpdir, 't');
system(sprintf('wb_command -cifti-separate %s_t_unthr.dscalar.nii COLUMN -volume-all %s.nii', figbasename, imgname));

cbvert = suit_map2surf(sprintf('%s.nii', imgname), 'space', 'SUIT', 'stats', @minORmax);
cbvert(cbvert > 0) = cbvert(cbvert > 0) ./ t_fdr_thr(4);
cbvert(cbvert < 0) = -cbvert(cbvert < 0) ./ t_fdr_thr(1);
tcol = [cols.fx_interp(cols.wb_rybcbl, 1000)];
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', tcol, 'cscale', [-1 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_t_unthr_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', cols.wb_power(ismember(parc.netnames, 'CB'), :), 'cscale', [-1 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_net_unthr_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

cbvert(cbvert < t_fdr_thr(3)/t_fdr_thr(4) & cbvert > -t_fdr_thr(2)/t_fdr_thr(1)) = NaN;
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', tcol, 'cscale', [-1 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_t_fdr05_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', cols.wb_power(ismember(parc.netnames, 'CB'), :), 'cscale', [-1 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_net_fdr05_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

rmdir(tmpdir, 's');

%% Visualize subcortex

tmpdir = tempname(tempdir); mkdir(tmpdir);

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.data(t_t < 0) = NaN;
tempcf.data(~ismember(parc.id, find(ismember(parc.net, find(ismember(parc.netnames, {'THA', 'BG'})))))) = NaN;
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 't_pos_thabg');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.data(t_t > 0) = NaN;
tempcf.data(~ismember(parc.id, find(ismember(parc.net, find(ismember(parc.netnames, {'THA', 'BG'})))))) = NaN;
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 't_neg_thabg');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.data(t_t < 0) = NaN;
tempcf.data(~ismember(parc.id, find(ismember(parc.net, find(ismember(parc.netnames, {'HCAMY', 'BS'})))))) = NaN;
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 't_pos_hcamybs');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = t_t;
tempcf.data(t_t > 0) = NaN;
tempcf.data(~ismember(parc.id, find(ismember(parc.net, find(ismember(parc.netnames, {'HCAMY', 'BS'})))))) = NaN;
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 't_neg_hcamybs');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

make_CLUT(cols.wb_rybcbl_rybl, fullfile(mricroglclutdir, 'wb_rybcbl_rybl.clut'));
make_CLUT(cols.wb_rybcbl_blbc, fullfile(mricroglclutdir, 'wb_rybcbl_blbc.clut'));

net_idx = find(ismember(parc.netnames, {'THA', 'HCAMY', 'BG', 'BS'})).';
for i = 1:numel(net_idx)
    tempcf = parc.cf;
    tempcf.data = ones(size(parc.id));
    tempcf.data(t_fdr_h) = 2;
    tempcf.data(~ismember(parc.id, find(parc.net == net_idx(i)))) = NaN;
    tempcf.mapname = {''};
    imgname = fullfile(tmpdir, sprintf('net%d', net_idx(i)));
    ft_write_cifti_mod(imgname, tempcf);
    system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));
    make_CLUT(cols.wb_power(net_idx(i), :), fullfile(mricroglclutdir, sprintf('wb_power_net%d.clut', net_idx(i))));
end

mricrogl_command = ['import gl', newline, ...
    'gl.resetdefaults()', newline, ...
    'gl.backcolor(255,255,255)', newline, ...
    'gl.colorbarposition(0)', newline, ...
    sprintf('gl.loadimage("%s/data/standard/MNI152_T1_1mm_brain.nii.gz")', getenv('FSLDIR')), newline, ...
    'gl.overlayloadsmooth(0)', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 't_pos_thabg')), newline, ...
    sprintf('gl.minmax(1, 0, %.4f)', t_fdr_thr(4)), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorfromzero(1, 1)', newline, ...
    'gl.colorname(1, "wb_rybcbl_rybl")', newline, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 't_neg_thabg')), newline, ...
    sprintf('gl.minmax(2, %.4f, 0)', t_fdr_thr(1)), newline, ...
    'gl.zerointensityinvisible(2, 1)', newline, ...
    'gl.colorfromzero(2, 1)', newline, ...
    'gl.colorname(2, "wb_rybcbl_blbc")', newline, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_t_unthr_z_18_9_-1.png")', figbasename), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', t_fdr_thr([3 4])), newline, ...
    sprintf('gl.minmax(2, %.4f, %.4f)', t_fdr_thr([1 2])), newline, ...
    'gl.zerointensityinvisible(1, 0)', newline, ...
    'gl.zerointensityinvisible(2, 0)', newline, ...
    'gl.invertcolor(2, 1)', newline, ...
    sprintf('gl.savebmp("%s_t_fdr05_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 't_pos_hcamybs')), newline, ...
    sprintf('gl.minmax(1, 0, %.4f)', t_fdr_thr(4)), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorfromzero(1, 1)', newline, ...
    'gl.colorname(1, "wb_rybcbl_rybl")', newline, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 't_neg_hcamybs')), newline, ...
    sprintf('gl.minmax(2, %.4f, 0)', t_fdr_thr(1)), newline, ...
    'gl.zerointensityinvisible(2, 1)', newline, ...
    'gl.colorfromzero(2, 1)', newline, ...
    'gl.colorname(2, "wb_rybcbl_blbc")', newline, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_t_unthr_x_-25_-2_25.png")', figbasename), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', t_fdr_thr([3 4])), newline, ...
    sprintf('gl.minmax(2, %.4f, %.4f)', t_fdr_thr([1 2])), newline, ...
    'gl.zerointensityinvisible(1, 0)', newline, ...
    'gl.zerointensityinvisible(2, 0)', newline, ...
    'gl.invertcolor(2, 1)', newline, ...
    sprintf('gl.savebmp("%s_t_fdr05_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

net_idx = find(ismember(parc.netnames, {'THA', 'BG'})).';
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, sprintf('net%d', net_idx(i)))), newline, ...
        sprintf('gl.minmax(%d, 0, 2)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 1)', i), newline, ...
        sprintf('gl.colorname(%d, "wb_power_net%d")', i, net_idx(i)), newline];
end
mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_net_unthr_z_18_9_-1.png")', figbasename), newline];
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.minmax(%d, 1, 2)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 0)', i), newline];
end
mricrogl_command = [mricrogl_command, ...
    sprintf('gl.savebmp("%s_net_fdr05_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

net_idx = find(ismember(parc.netnames, {'HCAMY', 'BS'})).';
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, sprintf('net%d', net_idx(i)))), newline, ...
        sprintf('gl.minmax(%d, 0, 2)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 1)', i), newline, ...
        sprintf('gl.colorname(%d, "wb_power_net%d")', i, net_idx(i)), newline];
end
mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_net_unthr_x_-25_-2_25.png")', figbasename), newline];
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.minmax(%d, 1, 2)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 0)', i), newline];
end
mricrogl_command = [mricrogl_command, ...
    sprintf('gl.savebmp("%s_net_fdr05_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, 'gl.quit()'];

mricrogl_py = fullfile(tmpdir, 'mricrogl.py');
fid = fopen(mricrogl_py, 'wt');
fwrite(fid, mricrogl_command);
fclose(fid);

system(sprintf('%s %s', mricrogldir, mricrogl_py));

rmdir(tmpdir, 's');

delete(fullfile(mricroglclutdir, 'wb_rybcbl_rybl.clut'));
delete(fullfile(mricroglclutdir, 'wb_rybcbl_blbc.clut'));

net_idx = find(ismember(parc.netnames, {'THA', 'HCAMY', 'BG', 'BS'})).';
for i = 1:numel(net_idx)
    delete(fullfile(mricroglclutdir, sprintf('wb_power_net%d.clut', net_idx(i))));
end

%% Prepare network

net_labels = {'DMN', 'AMTL', 'PMTL', 'CAN', 'MVIS', 'LVIS', ...
    'PMN', 'FPN', 'DAN', 'PMOT', 'LANG', 'SAL', 'AMN', ...
    'SCAN', 'FSMN', 'HSMN', 'LSMN', 'AUD', ...
    'CB', 'THA', 'BG', 'HCAMY', 'BS'}; % No NONE
net_labels = net_labels(ismember(net_labels, parc.netnames(vertnet(vertnet ~= 0))));
[~, net_vidx] = ismember(parc.netnames(vertnet(vertnet ~= 0)), net_labels);
net_c = cols.wb_power(vertnet(vertnet ~= 0), :);
net_c = net_c(net_vidx ~= 0, :);
net_x = net_vidx(net_vidx ~= 0);
net_y = t_t(vertnet ~= 0);
net_y = net_y(net_vidx ~= 0);
wh_thr = t_fdr_h(vertnet ~= 0);
wh_thr = wh_thr(net_vidx ~= 0);

%%

pie_dia = histcounts(net_x, 1:max(net_x)+1);
pie_dia = (pie_dia ./ max(pie_dia)) .^ 0.5;
for i = 1:max(net_x)
    colororder([0.9 0.9 0; 0.7 0.7 0.7; 0 0.9 0.9]);
    piechart([sum(net_y(net_x==i) >= t_fdr_thr(3)), ...
        sum(net_y(net_x==i) > t_fdr_thr(2) & net_y(net_x==i) < t_fdr_thr(3)), ...
        sum(net_y(net_x==i) <= t_fdr_thr(2))], 'Labels', {'', '', ''}, 'EdgeColor', 'none');
    
    resize_axes(gca, pie_dia(i) .* [800 800]);

    figname = sprintf('%s_t_net_pie_net%.2d.pdf', figbasename, i);
    exportgraphics(gcf, figname, 'BackgroundColor', 'none');
    close all;
end

%% Network-level swarm and box

hold on;
yregion(t_fdr_thr(2:3), 'FaceColor', [0.7 0.7 0.7]);
for i = 1:max(net_x)
    [f, xf] = kde(net_y(net_x==i));
    wh_real_xf = xf >= min(net_y(net_x==i)) & xf <= max(net_y(net_x==i));
    f = f(wh_real_xf);
    xf = xf(wh_real_xf);
    c = unique(net_c(net_x==i, :), 'rows');
    c = max(0, min(1, lab2rgb(rgb2lab(c) + [5 0 0; -5 0 0])));
    patch(i + [0;f;0;-flip(f)].*0.85, [xf(1);xf;xf(end);flip(xf)], c(1,:), 'LineWidth', 2, 'EdgeColor', c(2,:));
end
boxchart(net_x, net_y, 'MarkerStyle', 'none', 'LineWidth', 2, 'BoxWidth', 0.1, ...
    'BoxFaceAlpha', 1, 'BoxFaceColor', [1 1 1], 'BoxEdgeColor', [0.5 0.5 0.5], 'BoxMedianLineColor', [0.4 0.4 0.4], 'WhiskerLineColor', [0.6 0.6 0.6]);

set(gca, 'LineWidth', 2, 'XLim', [0 max(net_x)+1], 'XTick', 1:max(net_x), 'XTickLabel', net_labels, ...
    'TickLength', [0.005 0.005], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);
switch sj_num
    case 1
        set(gca, 'YLim', [-10 10], 'YTick', -10:2:10, 'YGrid', 'on', 'GridLineStyle', '--');
    case 2
        set(gca, 'YLim', [-8 8], 'YTick', -8:2:8, 'YGrid', 'on', 'GridLineStyle', '--');
end
set(get(gca, 'XAxis'), 'TickDirection', 'none');
set(get(gca, 'YAxis'), 'Exponent', false);

resize_axes(gca, [1341 250]);
set(gcf, 'color', 'w');

figname = sprintf('%s_t_net_violinbox.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname);
close all;

fprintf('Done.\n');

end
