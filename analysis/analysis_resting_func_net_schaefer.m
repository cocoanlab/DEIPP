function analysis_resting_func_net_schaefer(sj_num, bctgitdir, sctxparc, sctxparclabel, schaeferparc, schaeferlabel, varargin)

preptype = 'cocoan-preproc';
antype = 'cocoan-analysis';
bidsdir = fileparts(fileparts(mfilename('fullpath'))); % mfilename: bidsdir/code/~.m
if isempty(bidsdir); bidsdir = fileparts(pwd); end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'preptype'}
                preptype = varargin{i+1};
            case {'antype'}
                antype = varargin{i+1};
            case {'bidsdir'}
                bidsdir = varargin{i+1};
        end
    end
end

addpath(genpath(bctgitdir));
addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));
rcols = call_colors;

%% Basic setting

prepdir = fullfile(bidsdir, 'derivatives', preptype);
andir = fullfile(bidsdir, 'derivatives', antype);

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));

sj_id = tbl.participant_id(sj_num,:);
tasklabel = 'task-rest';

indnetdir = fullfile(andir, 'MSC_ind_net_parc', sj_id);
refnetdir = fullfile(andir, 'resting_func_net', sj_id);
if ~exist(refnetdir, 'dir'); mkdir(refnetdir); end

metadir = fullfile(prepdir, 'metadata');
anatmetafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, 'ses-01', 'anat'));
load(anatmetafile, 'anat');
anat = change_bidsdir(anat, bidsdir);
func = [];
funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabel));
funcmetafile = split(deblank(ls(funcmetafile)));
for ses_i = 1:numel(funcmetafile)
    tempstruct = load(funcmetafile{ses_i}, 'func');
    func{ses_i} = change_bidsdir(tempstruct.func, bidsdir);
end

phenodir = fullfile(bidsdir, 'phenotype');
sur_V1_tab = readtable(fullfile(phenodir, 'survey_V1.tsv'), 'FileType', 'text');
sur_V5_tab = readtable(fullfile(phenodir, 'survey_V5.tsv'), 'FileType', 'text');
sur_V10_tab = readtable(fullfile(phenodir, 'survey_V10.tsv'), 'FileType', 'text');

%% PC score of pain measurements

surtypes = {'VAS_D', 'VAS_WA', 'VAS_WW', 'SFMPQ_S', 'SFMPQ_A', 'PainDETECT', 'WPI', 'SSS'};
pain_mat = table2array(sur_V1_tab(strcmp(sur_V1_tab.participant_id, sj_id), surtypes));
[pain_pc,~,~,~,pain_pcexp] = pca(pain_mat);
pain_sc = pain_mat * pain_pc(:,1);
pain_sc = pain_sc(2:end);
pain_sc_norm = (pain_sc - min(pain_sc)) ./ (max(pain_sc) - min(pain_sc)) .* 100;

%% Prepare meta info of individual parcellation and ROI timeseries

schaefer.cf = ft_read_cifti_mod(func{1}.cf.bold);
schaefer.cf.data = [];
schaefer.cf.dimord = 'scalar_pos';

ctxparc_cf = ft_read_cifti_mod(schaeferparc);
ctxparc_cf.data(ctxparc_cf.data == 0) = NaN;
ctxparc_cf.data = ctxparc_cf.data(schaefer.cf.brainstructure(1:numel(ctxparc_cf.data)) ~= -1);

tmpdir = tempname(tempdir); mkdir(tmpdir);
system(sprintf('cp %s %s', func{1}.cf.bold, fullfile(tmpdir, 'template.dtseries.nii')));
system(sprintf('wb_command -cifti-create-dense-from-template %s %s -volume-all %s', ...
    fullfile(tmpdir, 'template.dtseries.nii'), fullfile(tmpdir, 'sctxatlas.dscalar.nii'), sctxparc));
sctxparc_cf = ft_read_cifti_mod(fullfile(tmpdir, 'sctxatlas.dscalar.nii'));
sctxparc_cf.data(sctxparc_cf.data == 0) = NaN;
sctxparc_cf.data = sctxparc_cf.data + max(ctxparc_cf.data);

brainstr_nomw = schaefer.cf.brainstructure(schaefer.cf.brainstructure ~= -1);
schaefer.id = zeros(numel(brainstr_nomw), 1);
schaefer.id(ismember(brainstr_nomw, 1:2)) = ctxparc_cf.data(ismember(brainstr_nomw, 1:2));
schaefer.id(~ismember(brainstr_nomw, 1:2)) = sctxparc_cf.data(~ismember(brainstr_nomw, 1:2));

schaefer.thrall = [0.003 0.004 0.005:0.005:0.050];
schaefer.thrref = 0.025;

s = load(sctxparclabel); f = fieldnames(s);
sctxparclabeldat = s.(f{1});

parcinfo = fileread(schaeferlabel);
parcinfo = reshape(split(deblank(parcinfo)), 6, []).';
parcnet = cellfun(@(a) a{3}, cellfun(@(a) split(a, '_'), parcinfo(:,1), 'un', false), 'un', false);

[~, ~, schaefer.net] = unique(parcnet, 'stable');
schaefer.net = [schaefer.net; repmat(sctxparclabeldat.dat(:,2)+max(schaefer.net), 1, size(schaefer.net,2))];

schaefer.netlabels = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN', 'THA', 'HCAMY', 'BG', 'CB', 'BS'};

fsLRdir = fullfile(prepdir, 'ciftify', anat.basename, 'MNINonLinear', 'fsaverage_LR32k');
surf_L = fullfile(fsLRdir, sprintf('%s.L.midthickness.32k_fs_LR.surf.gii', anat.basename));
surf_R = fullfile(fsLRdir, sprintf('%s.R.midthickness.32k_fs_LR.surf.gii', anat.basename));
schaefer.cf.pos(1:64984,:) = [getfield(gifti(surf_L), 'vertices'); getfield(gifti(surf_R), 'vertices')];

aparcfile = fullfile(fsLRdir, sprintf('%s.aparc.32k_fs_LR.dlabel.nii', anat.basename));
aparccf = ft_read_cifti_mod(aparcfile);
schaefer.novis = ~ismember(schaefer.net, 1) & ...
    ~ismember(1:numel(schaefer.net), unique(schaefer.id(find(ismember(aparccf.data, [5 7 11 13 21 40 42 46 48 56]))))).';

save(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_meta.mat', sj_id)), 'schaefer');

tempcf = schaefer.cf;
tempcf.data = schaefer.id;
tempcf.mapname = {'ROI number'};
imgname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_roi', sj_id));
ft_write_cifti_mod(imgname, tempcf);

tempcf = schaefer.cf;
tempcf.data = schaefer.id .* 0 + 1;
tempcf.mapname = {'ROI existence'};
imgname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_binaryroi', sj_id));
ft_write_cifti_mod(imgname, tempcf);

tempcf = schaefer.cf;
tempcf.data = schaefer.id;
for i = 1:max(schaefer.id)
    tempcf.data(schaefer.id == i) = schaefer.net(i);
end
tempcf.mapname = {'Network labels'};
imgname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_net', sj_id));
ft_write_cifti_mod(imgname, tempcf);

tempcf = schaefer.cf;
tempcf.data = schaefer.cf.brainstructure(schaefer.cf.brainstructure ~= -1);
tempcf.mapname = {'Brainstructure labels'};
imgname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_brainstructure', sj_id));
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = schaefer.cf;
tempcf.data = double(ismember(schaefer.id, find(~schaefer.novis)));
tempcf.mapname = {'ROI in the visual network and occipital areas'};
imgname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_binaryroivis', sj_id));
ft_write_cifti_mod(imgname, tempcf);

rf_ts = [];
for ses_i = 1:numel(funcmetafile)    
    tempcf = ft_read_cifti_mod(func{ses_i}.cf.bold);
    rf_dat = double(tempcf.data(:, logical(importdata(func{ses_i}.den.cen)))).';
    rf_dat = rf_dat - mean(rf_dat);
    rf_ts{ses_i} = splitapply(@(x) mean(x,2), rf_dat, schaefer.id.');
end
save(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_timeseries.mat', sj_id)), 'rf_ts');

%% Prepare connectivity

load(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_meta.mat', sj_id)), 'schaefer');
load(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_timeseries.mat', sj_id)), 'rf_ts');

rf_corr = cellfun(@(x) atanh(corrcoef(x)), rf_ts, 'un', false);
rf_corr_all = atanh(corrcoef(cat(1, rf_ts{:})));

assert(~any(cellfun(@(x) any(isnan(x), 1:2), cat(2, rf_corr, rf_corr_all))));

%% session-level similarity

r = cellfun(@refmt_r, rf_corr, 'un', false);
r = corrcoef(cat(2, r{:}));
imagesc(r, [0.45 0.85]);
colormap(rcols.bipolar);
colorbar('TickLength', 0);
seslabels = strcat('ses-', split(deblank(sprintf('%.2d\n', 1:size(r,1)))));
set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out', ...
    'XTick', 1:size(r,1), 'XTickLabel', seslabels, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(r,1), 'YTickLabel', seslabels, 'YTickLabelRotation', 0);
resize_axes(gca, repmat(size(r,1) * 25, 1, 2));
set(gcf, 'Color', 'w');
figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_corr_session_similarity.pdf', sj_id));
pagesetup(gcf); saveas(gcf, figname);
close all;

%% Network attributes

Eglob = NaN(numel(schaefer.thrall), numel(rf_corr));
ModQ = NaN(numel(schaefer.thrall), numel(rf_corr));
for thr_i = 1:numel(schaefer.thrall)
    for ses_i = 1:numel(rf_corr)
        r = threshold_proportional(tanh(rf_corr{ses_i}), schaefer.thrall(thr_i));
        Eglob(thr_i,ses_i) = efficiency_wei(r);
        ModQ(thr_i,ses_i) = calc_modularity(r, schaefer.net);
    end
end
save(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_netattr.mat', sj_id)), 'Eglob', 'ModQ');

%% Correlation: Network attribute - Pain PC

load(fullfile(refnetdir, sprintf('%s_resting_func_schaefer_netattr.mat', sj_id)), 'Eglob', 'ModQ');

for netattr = {Eglob ModQ; 'Eglob' 'ModQ'; 'Global efficiency' 'Modularity Q'}
    [x, y] = deal(mean(netattr{1}).', pain_sc);
    [r, p] = corr(x, y);
    fprintf('Correlation: r = %.2f, p = %.4f\n', r, p);
    scatter(x, y, 50, [0.6 0.6 0.6], 'filled');
    xlabel(netattr{3});
    ylabel('Normalized Pain PC score');
    set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out', ...
        'XLim', extlim(x, 0.1), 'YLim', extlim(y, 0.1));
    refline;
    resize_axes(gca, [200 200]);
    set(gcf, 'Color', 'w');
    figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_corr_%s_painPC.pdf', sj_id, netattr{2}));
    pagesetup(gcf); saveas(gcf, figname);
    close all;
end

%% Connectivity matrix

r = tanh(rf_corr_all);
vis_grp_r(r, schaefer.net, rcols.power);
clim_val = extlim(r);
fprintf('Clim: [%.4f, %.4f]\n', clim_val(:));
caxis(clim_val);
colormap(rcols.bipolar);
colorbar('TickLength', 0);
resize_axes(gca, [700, 700]);
figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_rawconn_ses-all.pdf', sj_id));
pagesetup(gcf); saveas(gcf, figname);
close all;

r = tanh(grpmean_r(rf_corr_all, schaefer.net));
vis_grp_r(r, [], rcols.power(unique(schaefer.net),:));
clim_val = extlim(r);
fprintf('Clim: [%.4f, %.4f]\n', clim_val(:));
caxis(clim_val);
colormap(rcols.bipolar);
colorbar('TickLength', 0);
resize_axes(gca, [500, 500]);
figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_netconn_ses-all.pdf', sj_id));
pagesetup(gcf); saveas(gcf, figname);
close all;

%% Correlation: connectivity - Pain PC

x = cellfun(@refmt_r, rf_corr, 'un', false);
x = cat(2, x{:}).';
y = pain_sc;
[r, p] = corr(x, y);
r = refmt_r(r);
vis_grp_r(r, schaefer.net, rcols.power);
clim_val = [-1 1] .* abs(max(r(:)));
fprintf('Clim: [%.4f, %.4f]\n', clim_val(:));
caxis(clim_val);
colormap(rcols.bipolar);
colorbar('TickLength', 0);
resize_axes(gca, [700, 700]);
figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_corr_rawconn_painPC.pdf', sj_id));
pagesetup(gcf); saveas(gcf, figname);
close all;

x = cellfun(@(x) grpmean_r(x, schaefer.net), rf_corr, 'un', false);
x = cellfun(@refmt_r, x, 'un', false);
x = cat(2, x{:}).';
y = pain_sc;
[r, p] = corr(x, y);
r = refmt_r(r);
vis_grp_r(r, [], rcols.power(unique(schaefer.net),:));
sig_unc = refmt_r(p < 0.05);
[sig_unc_row, sig_unc_col] = find(sig_unc);
for sig_i = 1:numel(sig_unc_row)
    rectangle('Position', [sig_unc_row(sig_i)-0.45 sig_unc_col(sig_i)-0.45 0.9 0.9], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 4);
end
sig_fdr = refmt_r(p <= max([FDR(p,0.05), 0]));
[sig_fdr_row, sig_fdr_col] = find(sig_fdr);
for sig_i = 1:numel(sig_fdr_row)
    rectangle('Position', [sig_fdr_row(sig_i)-0.45 sig_fdr_col(sig_i)-0.45 0.9 0.9], 'EdgeColor', [0.95 0.95 0.95], 'LineWidth', 4);
end
clim_val = [-1 1] .* abs(max(r(:)));
fprintf('Clim: [%.4f, %.4f]\n', clim_val(:));
caxis(clim_val);
colormap(rcols.bipolar);
colorbar('TickLength', 0);
resize_axes(gca, [500, 500]);
figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_corr_netconn_painPC.pdf', sj_id));
pagesetup(gcf); saveas(gcf, figname);
close all;

%% spring network

% NONE DMN LVN FPN MVN DAN PM VAN SAN CON HSMN FSMN AN AMTL PMTL PAM CT LSMN THA HCAMY BG CB BS
% DMN/AMTL/PMTL - FPN/PAM - DAN - SAN/CON - VAN/CT
% 135 67.5 0 270 202.5

% xy_net = deg2rad([150 30 90 30 60 270 180 210 210 270 270 240 150 150 90 180 270 NaN].');
xy_net = deg2rad([NaN 135 NaN 67.5 NaN 0 NaN 202.5 270 270 NaN NaN NaN 135 135 67.5 202.5 NaN NaN NaN NaN NaN NaN].');
xy_net = [cos(xy_net) sin(xy_net)] .* 2 + [-1 1];
xy_net([3 5 7 11 12 13 18], :) = repmat([1 -1], 7, 1);
xy_net([1 19 20 21 22 23], :) = repmat([0 0], 6, 1);
xy_ref = repmat({xy_net(schaefer.net, 1:2)}, 1, numel(schaefer.thrall)+1);

for thr_i = 1:numel(schaefer.thrall)
    A = threshold_proportional(tanh(rf_corr_all), schaefer.thrall(thr_i));
    [A_comps, A_comp_sizes] = get_components(A);
    wh_largecomp = ismember(A_comps, find(A_comp_sizes > floor(size(A,1) * 0.05)));
    G = graph(A(wh_largecomp, wh_largecomp));
    h_graph = plot(G, 'Layout', 'force', ...
        'WeightEffect', 'inverse', 'UseGravity', true, 'Iterations', 1000, ...
        'XStart', xy_ref{thr_i}(wh_largecomp,1), 'YStart', xy_ref{thr_i}(wh_largecomp,2));
    xy_ref{thr_i+1}(wh_largecomp, :) = [h_graph.XData.' h_graph.YData.'];
    h_graph_xy = {h_graph.XData, h_graph.YData};
    h_graph_xylim = {xlim, ylim};
    close all; figure; hold on;
    line(h_graph_xy{1}(G.Edges.EndNodes).', h_graph_xy{2}(G.Edges.EndNodes).', 'Color', [0.6 0.6 0.6 0.7]);
    scatter(h_graph_xy{1}, h_graph_xy{2}, 100, rcols.power(schaefer.net(wh_largecomp), :), 'MarkerFaceColor', 'flat', ...
        'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1);
    set(gca, 'XLim', h_graph_xylim{1}, 'YLim', h_graph_xylim{2}, 'Visible', 'off');
    set(gcf, 'Position', [420         104        1077         875], 'Color', 'w');
    annotation(gcf, 'textbox', [0 0 1 1], 'String', sprintf('%s\nTop %.1f%%', sj_id, schaefer.thrall(thr_i)*100), 'EdgeColor', 'none', 'FontSize', 16);
    
    pagesetup(gcf);
    figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_spring_embed_%s.pdf', sj_id, strrep(sprintf('%.3f', schaefer.thrall(thr_i)),'.','_')));
    saveas(gcf, figname);
    figname = fullfile(refnetdir, sprintf('%s_resting_func_schaefer_spring_embed_%s.png', sj_id, strrep(sprintf('%.3f', schaefer.thrall(thr_i)),'.','_')));
    saveas(gcf, figname);
    close all;
end

for ses_i = 1:numel(rf_corr)
    A = threshold_proportional(tanh(rf_corr{ses_i}), schaefer.thrref);
    [A_comps, A_comp_sizes] = get_components(A);
    wh_largecomp = ismember(A_comps, find(A_comp_sizes > floor(size(A,1) * 0.05)));
    G = graph(A(wh_largecomp, wh_largecomp));
    h_graph = plot(G, 'Layout', 'force', ...
        'WeightEffect', 'inverse', 'UseGravity', true, 'Iterations', 1000, ...
        'XStart', xy_ref{find(schaefer.thrall==schaefer.thrref)+1}(wh_largecomp,1), 'YStart', xy_ref{find(schaefer.thrall==schaefer.thrref)+1}(wh_largecomp,2));
    h_graph_xy = {h_graph.XData, h_graph.YData};
    h_graph_xylim = {xlim, ylim};
    close all; figure; hold on;
    line(h_graph_xy{1}(G.Edges.EndNodes).', h_graph_xy{2}(G.Edges.EndNodes).', 'Color', [0.6 0.6 0.6 0.7]);
    scatter(h_graph_xy{1}, h_graph_xy{2}, 100, rcols.power(schaefer.net(wh_largecomp), :), 'MarkerFaceColor', 'flat', ...
        'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1);
    set(gca, 'XLim', h_graph_xylim{1}, 'YLim', h_graph_xylim{2}, 'Visible', 'off');
    set(gcf, 'Position', [420         104        1077         875], 'Color', 'w');
    annotation(gcf, 'textbox', [0 0 1 1], 'String', sprintf('%s, ses-%.2d\nTop %.1f%%', sj_id, ses_i, schaefer.thrref*100), 'EdgeColor', 'none', 'FontSize', 16);
    
    fig_pos = get(gcf, 'Position');
    
    h_tab_netattr = uitable('Data', [mean(Eglob(:,ses_i)), mean(ModQ(:,ses_i))], ...
        'Rowname', '', 'ColumnName', {'Eglob', 'ModQ'}, 'FontSize', 14, 'Position', [0 0 1000 300]);
    h_tab_netattr.Position(3:4) = h_tab_netattr.Extent(3:4);
    h_tab_netattr.Position(1) = fig_pos(3) - (h_tab_netattr.Position(1) + h_tab_netattr.Position(3));
    
    wh_sjses_V1 = strcmp(sur_V1_tab.participant_id, sj_id) & strcmp(sur_V1_tab.session_id, sprintf('ses-%.2d', ses_i));
    varname_V1 = {'VAS_D', 'VAS_WA', 'VAS_WW', 'SFMPQ_S', 'SFMPQ_A', 'PainDETECT', 'WPI', 'SSS', 'HADS_A', 'HADS_D'};
    h_tab_V1 = uitable('Data', [round(pain_sc_norm(ses_i)) table2array(sur_V1_tab(wh_sjses_V1, varname_V1))], ...
        'Rowname', '', 'ColumnName', [{'Pain PC'}, varname_V1], 'FontSize', 14, 'Position', [0 0 1000 300]);
    h_tab_V1.Position(3:4) = h_tab_V1.Extent(3:4);
    h_tab_V1.Position(1:2) = fig_pos(3:4) - (h_tab_V1.Position(1:2) + h_tab_V1.Position(3:4));
    if ~mod(ses_i, 5)
        wh_sjses_V5 = strcmp(sur_V5_tab.participant_id, sj_id) & strcmp(sur_V5_tab.session_id, sprintf('ses-%.2d', ses_i));
        varname_V5 = {'FIQ', 'PDI', 'TSK', 'PANAS_P', 'PANAS_N'};
        h_tab_V5 = uitable('Data', table2array(sur_V5_tab(wh_sjses_V5, varname_V5)), ...
            'Rowname', '', 'ColumnName', varname_V5, 'FontSize', 14, 'Position', [0 h_tab_V1.Position(4) 1000 300]);
        h_tab_V5.Position(3:4) = h_tab_V5.Extent(3:4);
        h_tab_V5.Position(1:2) = h_tab_V1.Position(1:2) - [h_tab_V5.Position(3)-h_tab_V1.Position(3) h_tab_V5.Position(4)];
    end
    if ~mod(ses_i, 10)
        wh_sjses_V10 = strcmp(sur_V10_tab.participant_id, sj_id) & strcmp(sur_V10_tab.session_id, sprintf('ses-%.2d', ses_i));
        varname_V10 = {'PCS_T', 'PCS_R', 'PCS_M', 'PCS_H', 'CPVI_I', 'CPVI_S', 'CPVI_D', 'CFQ', 'CPAQ_T', 'CPAQ_AE', 'CPAQ_PW'};
        h_tab_V10 = uitable('Data', table2array(sur_V10_tab(wh_sjses_V10, varname_V10)), ...
            'Rowname', '', 'ColumnName', varname_V10, 'FontSize', 14, 'Position', [0 h_tab_V5.Position(4) 1000 300]);
        h_tab_V10.Position(3:4) = h_tab_V10.Extent(3:4);
        h_tab_V10.Position(1:2) = h_tab_V5.Position(1:2) - [h_tab_V10.Position(3)-h_tab_V5.Position(3) h_tab_V10.Position(4)];
    end
    
    pagesetup(gcf);
    figname = fullfile(refnetdir, sprintf('%s_ses-%.2d_resting_func_schaefer_spring_embed_%s.pdf', sj_id, ses_i, strrep(sprintf('%.3f', schaefer.thrref),'.','_')));
    saveas(gcf, figname);
    figname = fullfile(refnetdir, sprintf('%s_ses-%.2d_resting_func_schaefer_spring_embed_%s.png', sj_id, ses_i, strrep(sprintf('%.3f', schaefer.thrref),'.','_')));
    saveas(gcf, figname);
    close all;
end

end

