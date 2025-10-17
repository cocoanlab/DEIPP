function analysis_visualize_model_netimpwei(sj_num, spmdir, mricrogldir, varargin)

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
preddir = fullfile(andir, 'predict_rating', sj_id);

cols = call_colors;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

featimpname = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s.mat', sj_id, parctype));
if ~exist(featimpname, 'file'); return; end
load(featimpname, 'fimp');

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

figbasename = fullfile(preddir, sprintf('%s_predict_rating_netimpwei_%s', sj_id, parctype));

%% Prepare network, feature importance, and model weight

w = refmt_r(out.beta_best_edge.all(2:end));
wpos = sum(w .* double(w > 0), 2);
wneg = sum(w .* double(w < 0), 2);

net_labels = {'DMN', 'AMTL', 'PMTL', 'CAN', ...
    'PMN', 'FPN', 'DAN', 'PMOT', 'LANG', 'SAL', 'AMN', ...
    'SCAN', 'FSMN', 'HSMN', 'LSMN', 'AUD', ...
    'CB', 'THA', 'BG', 'HCAMY', 'BS'}; % No NONE
net_labels = net_labels(ismember(net_labels, parc.netnames(parc.net(parc.whincl))));
[~, net_ridx] = ismember(parc.netnames(parc.net(parc.whincl)), net_labels);
net_c = cols.wb_power(parc.net(parc.whincl), :);
net_c = net_c(net_ridx ~= 0, :);
net_x = net_ridx(net_ridx ~= 0);
net_y = struct('imp', fimp.reg(net_ridx ~= 0), 'pos', wpos(net_ridx ~= 0), 'neg', wneg(net_ridx ~= 0));
topk = 10;
wh_topk = struct('imp', net_y.imp >= quantile(net_y.imp, 1 - topk/numel(net_y.imp)), ...
    'pos', net_y.pos >= quantile(net_y.pos, 1 - topk/numel(net_y.pos)), ...
    'neg', net_y.neg < quantile(net_y.neg, topk/numel(net_y.neg)));
topk_id = find(parc.whincl);
topk_id = topk_id(net_ridx ~= 0);
topk_id = topk_id(any([wh_topk.imp, wh_topk.pos, wh_topk.neg], 2));

%% Network-level swarm and box

for f = fieldnames(net_y).'
    hold on;
    boxchart(net_x, net_y.(f{1}), 'MarkerStyle', 'none', 'LineWidth', 2, 'BoxWidth', 0.4, ...
        'BoxFaceColor', 'none', 'BoxEdgeColor', [0.5 0.5 0.5], 'BoxMedianLineColor', [0.4 0.4 0.4], 'WhiskerLineColor', [0.6 0.6 0.6]);
    swarmchart(net_x, net_y.(f{1}), 60, net_c, 'filled', ...
        'AlphaData', double(wh_topk.(f{1})) .* 0.95, 'MarkerFaceAlpha', 'flat', 'MarkerEdgeAlpha', 'flat', ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.5, 'XJitterWidth', 2.4);
    swarmchart(net_x, net_y.(f{1}), 30, net_c, 'filled', ...
        'AlphaData', double(~wh_topk.(f{1})) .* 0.75, 'MarkerFaceAlpha', 'flat', 'MarkerEdgeAlpha', 'flat', ...
        'MarkerEdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.5, 'XJitterWidth', 0.8);
    
    set(gca, 'LineWidth', 2, 'XLim', [0 max(net_x)+1], 'XTick', 1:max(net_x), 'XTickLabel', net_labels, ...
        'TickLength', [0.005 0.005], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);
    switch sj_num
        case 1
            switch f{1}
                case 'imp'
                    set(gca, 'YLim', [-0.0336 0.05], 'YTick', -0.03:0.01:0.05, 'YGrid', 'on', 'GridLineStyle', '--');
                case 'pos'
                    set(gca, 'YLim', [0 0.08], 'YTick', 0:0.01:0.08, 'YGrid', 'on', 'GridLineStyle', '--');
                case 'neg'
                    set(gca, 'YLim', [-0.07 0], 'YTick', -0.07:0.01:0, 'YGrid', 'on', 'GridLineStyle', '--');
            end
        case 2
            switch f{1}
                case 'imp'
                    set(gca, 'YLim', [-0.0066 0.025], 'YTick', -0.005:0.005:0.025, 'YGrid', 'on', 'GridLineStyle', '--');
                case 'pos'
                    set(gca, 'YLim', [0 0.12], 'YTick', 0:0.02:0.12, 'YGrid', 'on', 'GridLineStyle', '--');
                case 'neg'
                    set(gca, 'YLim', [-0.16 0], 'YTick', -0.16:0.02:0, 'YGrid', 'on', 'GridLineStyle', '--');
            end
    end
    if strcmp(f{1}, 'neg'); set(gca, 'YDir', 'reverse'); end
    set(get(gca, 'XAxis'), 'TickDirection', 'none');
    set(get(gca, 'YAxis'), 'Exponent', false);
    
    resize_axes(gca, [1310 250]);
    set(gcf, 'color', 'w');
    
    figname = sprintf('%s_%s_swarmbox.pdf', figbasename, f{1});
    pagesetup(gcf); saveas(gcf, figname);
    close all;
end

%% Visualize cortex

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
tempcf.data(~isnan(parc.id)) = parc.netorig(parc.id(~isnan(parc.id)));
tempcf.data(~ismember(parc.id, topk_id)) = NaN;
tempcf.mapname = {'ROI network label'};
imgname = sprintf('%s_top%d', figbasename, topk);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));
% Participant 1: 55 -45 -135 0.7
% Participant 2: -45 -100 -100 0.7

%% Visualize cerebellum

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
cbvert = NaN(size(cbparc.cdata));
roi_idx = 1:max(cbparc.cdata);
for j = 1:numel(roi_idx)
    cbvert(cbparc.cdata == roi_idx(j)) = 1 - (j-1)/numel(roi_idx);
end
cbvert(~ismember(cbparc.cdata, find(ismember(find(strcmp(parc.netnames(parc.net), 'CB')), topk_id)))) = NaN;
cbvert(cbvert == 0) = NaN;
netcol = lab2rgb(rgb2lab(cols.wb_power(ismember(parc.netnames, 'CB'), :)) + [linspace(-3, 3, 255).' zeros(255, 2)]);
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', netcol, 'cscale', [0 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_top%d_cerebellum.png', figbasename, topk);
exportgraphics(gcf, figname);
close all;

%% Visualize subcortex

tmpdir = tempname(tempdir); mkdir(tmpdir);

net_idx = find(ismember(parc.netnames, {'THA', 'HCAMY', 'BG', 'BS'})).';
for i = 1:numel(net_idx)
    tempcf = parc.cf;
    tempcf.data = NaN(size(parc.id));
    tempcf.mapname = {''};
    roi_idx = find(parc.net == net_idx(i)).';
    nh_names = regexprep(parc.names(roi_idx), '-[lr]h$', '');
    [~, u_nh_idx] = unique(nh_names, 'stable');
    r_nh_idx = setdiff(1:numel(roi_idx), u_nh_idx);
    [~, r2u_nh_idx] = ismember(nh_names(r_nh_idx), nh_names(u_nh_idx));
    for j = 1:numel(u_nh_idx)
        tempcf.data(parc.id == roi_idx(u_nh_idx(j))) = 1 - (j-1)/numel(u_nh_idx);
    end
    for j = 1:numel(r_nh_idx)
        tempcf.data(parc.id == roi_idx(r_nh_idx(j))) = 1 - (r2u_nh_idx(j)-1)/numel(u_nh_idx) - 0.08;
    end
    tempcf.data(~ismember(parc.id, topk_id)) = NaN;
    imgname = fullfile(tmpdir, sprintf('roi_net%d', net_idx(i)));
    ft_write_cifti_mod(imgname, tempcf);
    system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));
    netcol = lab2rgb(rgb2lab(cols.wb_power(net_idx(i), :)) + [linspace(-3, 3, 255).' zeros(255, 2)]);
    make_CLUT(max(0, min(1, netcol)), fullfile(mricroglclutdir, sprintf('wb_power_net%d.clut', net_idx(i))));
end

mricrogl_command = ['import gl', newline, ...
    'gl.resetdefaults()', newline, ...
    'gl.backcolor(255,255,255)', newline, ...
    'gl.colorbarposition(0)', newline, ...
    sprintf('gl.loadimage("%s/data/standard/MNI152_T1_1mm_brain.nii.gz")', getenv('FSLDIR')), newline, ...
    'gl.overlayloadsmooth(0)', newline];

net_idx = find(ismember(parc.netnames, {'THA', 'BG'})).';
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, sprintf('roi_net%d', net_idx(i)))), newline, ...
        sprintf('gl.minmax(%d, 0, 1)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 1)', i), newline, ...
        sprintf('gl.colorname(%d, "wb_power_net%d")', i, net_idx(i)), newline];
end
mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_top%d_z_18_9_-1.png")', figbasename, topk), newline, ...
    'gl.overlaycloseall()', newline];

net_idx = find(ismember(parc.netnames, {'HCAMY', 'BS'})).';
for i = 1:numel(net_idx)
    mricrogl_command = [mricrogl_command, ...
        sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, sprintf('roi_net%d', net_idx(i)))), newline, ...
        sprintf('gl.minmax(%d, 0, 1)', i), newline, ...
        sprintf('gl.zerointensityinvisible(%d, 1)', i), newline, ...
        sprintf('gl.colorname(%d, "wb_power_net%d")', i, net_idx(i)), newline];
end
mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_top%d_x_-25_-2_25.png")', figbasename, topk), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, 'gl.quit()'];

mricrogl_py = fullfile(tmpdir, 'mricrogl.py');
fid = fopen(mricrogl_py, 'wt');
fwrite(fid, mricrogl_command);
fclose(fid);

system(sprintf('%s %s', mricrogldir, mricrogl_py));

rmdir(tmpdir, 's');

net_idx = find(ismember(parc.netnames, {'THA', 'HCAMY', 'BG', 'BS'})).';
for i = 1:numel(net_idx)
    delete(fullfile(mricroglclutdir, sprintf('wb_power_net%d.clut', net_idx(i))));
end

fprintf('Done.\n');

end
