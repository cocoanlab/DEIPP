function analysis_visualize_model_weights(sj_num, opt, spmdir, varargin)

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

addpath(genpath(spmdir));
addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
preddir = fullfile(andir, 'predict_rating_func', sj_id);
refnetdir = fullfile(andir, 'resting_func_net', sj_id);

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

anatmetafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, 'ses-01', 'anat'));
load(anatmetafile, 'anat');
anat = change_bidsdir(anat, bidsdir);
fsLRdir = fullfile(prepdir, 'ciftify', anat.basename, 'MNINonLinear', 'fsaverage_LR32k');

rcols = call_colors;

%% Load data

parcname = fullfile(refnetdir, sprintf('%s_resting_func_%s_meta.mat', sj_id, opt.parc));
if ~exist(parcname, 'file'); return; end
parc = load(parcname, opt.parc);
parc = parc.(opt.parc);

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s_%s.mat', ...
    sj_id, opt.dfc, opt.parc));
if ~exist(featname, 'file'); return; end
load(featname, 'roi_id', 'roi_net', 'wh_exclude');

featimpname = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
if ~exist(featimpname, 'file'); return; end
load(featimpname, 'fimp');

figbasename = strrep(strrep(featimpname, '_predict_rating_featimp_', '_predict_rating_figure_'), '.mat', '');

%% Feature importance

add_nonzero = eps;
imp_corr = fimp.reg.corr.mean;
imp_corr_scaled = imp_corr .* 1000;
imp_corr_scaled_minmax = [add_nonzero round(max(imp_corr_scaled), -1)];
imp_corr_scaled = min(max(imp_corr_scaled, imp_corr_scaled_minmax(1)), imp_corr_scaled_minmax(2));
fprintf('Original scale: %.4f - %.4f\n', min(imp_corr), max(imp_corr));
fprintf('Rearranged scale: %.4f - %.4f\n', imp_corr_scaled_minmax);

idx_exclude = find(wh_exclude);
[imp_corr_top, wh_imp_corr_top] = maxk(imp_corr, 5);
idx_imp_corr_top = idx_exclude(wh_imp_corr_top);
fprintf('Top 5 feature importance:\n\n');
disp(array2table([idx_imp_corr_top, imp_corr_top], 'VariableNames', {'ROI', 'Imp'}));

barcols = interp1(linspace(imp_corr_scaled_minmax(1),imp_corr_scaled_minmax(2),size(rcols.wb_cw_w,1)), ...
    rcols.wb_cw_w, imp_corr_scaled(wh_imp_corr_top));
h = barh(imp_corr_top, 0.3, 'EdgeColor', 'flat', 'FaceColor', 'flat', 'CData', barcols);
set(gca, 'YDir', 'reverse', 'LineWidth', 2, 'YLim', [0.65 numel(imp_corr_top)+0.35], 'YTick', [], ...
    'TickLength', [0.015 0.015], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);
switch sj_num
    case 1
        set(gca, 'XTick', 0:0.02:0.04);
    case 2
        set(gca, 'XTick', 0:0.01:0.02);
end
set(gcf, 'color', 'w');
resize_axes(gca, [130 398]);
figname = sprintf('%s_featimpcorr_top_bar.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname); saveas(gcf, strrep(figname, '.pdf', '.png'));
close all;

%% generate CIFTI

impcorrcf = parc.cf;
impcorrcf.data = NaN(size(roi_id));
for i = 1:numel(idx_exclude)
    impcorrcf.data(roi_id == idx_exclude(i)) = imp_corr(i);
end
impcorrcf.data = [impcorrcf.data NaN(size(impcorrcf.data,1), numel(idx_imp_corr_top))];
for i = 1:numel(idx_imp_corr_top)
    impcorrcf.data(roi_id == idx_imp_corr_top(i), i+1) = impcorrcf.data(roi_id == idx_imp_corr_top(i), 1);
end
impcorrcf.mapname = [{'Importance: Corr'}, ...
    strcat({'Importance: Corr, Top '}, num2str([1:numel(idx_imp_corr_top)].'), {' / '}, num2str(numel(idx_exclude))).'];
imgname = sprintf('%s_featimpcorr', figbasename);
ft_write_cifti_mod(imgname, impcorrcf);

impcorrcf_scaled = impcorrcf;
impcorrcf_scaled.data = NaN(size(roi_id));
for i = 1:numel(idx_exclude)
    impcorrcf_scaled.data(roi_id == idx_exclude(i)) = imp_corr_scaled(i);
end
impcorrcf_scaled.data = [impcorrcf_scaled.data NaN(size(impcorrcf_scaled.data,1), numel(idx_imp_corr_top))];
for i = 1:numel(idx_imp_corr_top)
    impcorrcf_scaled.data(roi_id == idx_imp_corr_top(i), i+1) = impcorrcf_scaled.data(roi_id == idx_imp_corr_top(i), 1);
end
impcorrcf_scaled.mapname = [{'Importance: Corr, scaled'}, ...
    strcat({'Importance: Corr, scaled, Top '}, num2str([1:numel(idx_imp_corr_top)].'), {' / '}, num2str(numel(idx_exclude))).'];
imgname = sprintf('%s_featimpcorr_scaled', figbasename);
ft_write_cifti_mod(imgname, impcorrcf_scaled);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

impcorrcf_thalbg = impcorrcf_scaled;
wh_thalbg = find(ismember(parc.net, find(contains(parc.netlabels, {'THA', 'BG'}))));
impcorrcf_thalbg.data = impcorrcf_scaled.data .* double(ismember(parc.id, wh_thalbg));
imgname = sprintf('%s_featimpcorr_scaled_thalbg', figbasename);
ft_write_cifti_mod(imgname, impcorrcf_thalbg);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

impcorrcf_hcamybs = impcorrcf_scaled;
wh_hcamybs = find(ismember(parc.net, find(contains(parc.netlabels, {'HCAMY', 'BS'}))));
impcorrcf_hcamybs.data = impcorrcf_scaled.data .* double(ismember(parc.id, wh_hcamybs));
imgname = sprintf('%s_featimpcorr_scaled_hcamybs', figbasename);
ft_write_cifti_mod(imgname, impcorrcf_hcamybs);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

%% Visualize cerebellum

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
wh_cb = roi_net == find(strcmp(parc.netlabels, 'CB'));
wh_cb = find(wh_cb(wh_exclude));
imp_corr_scaled_cb = imp_corr_scaled(wh_cb);

cbvert = NaN(size(cbparc.cdata));
for i = 1:max(cbparc.cdata)
    cbvert(cbparc.cdata == i) = 0.8 + i/max(cbparc.cdata)*0.4;
end
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', gray, 'cscale', [0 2], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_%s_cerebellum.png', figbasename, opt.parc);
exportgraphics(gcf, figname);
close all;

cbvert = NaN(size(cbparc.cdata));
for i = 1:max(cbparc.cdata)
    cbvert(cbparc.cdata == i) = imp_corr_scaled_cb(i);
end
cbvert(cbvert == 0) = NaN;
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', rcols.wb_cw_w, 'cscale', imp_corr_scaled_minmax, 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_featimpcorr_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

%% Visualize subcortex

mricrogl_command = ['import gl', newline, ...
    'gl.resetdefaults()', newline, ...
    'gl.backcolor(255,255,255)', newline, ...
    'gl.colorbarposition(0)', newline, ...
    sprintf('gl.loadimage("%s/data/standard/MNI152_T1_1mm_brain.nii.gz")', getenv('FSLDIR')), newline, ...
    'gl.overlayloadsmooth(0)', newline];

mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    ...
    sprintf('gl.overlayload("%s_featimpcorr_scaled_hcamybs.nii.gz")', figbasename), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', imp_corr_scaled_minmax), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_cw_w")', newline, ...
    sprintf('gl.savebmp("%s_featimpcorr_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    ...
    sprintf('gl.overlayload("%s_featimpcorr_scaled_thalbg.nii.gz")', figbasename), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', imp_corr_scaled_minmax), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_cw_w")', newline, ...
    sprintf('gl.savebmp("%s_featimpcorr_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, 'gl.quit()'];

tmpdir = tempname(tempdir); mkdir(tmpdir);
mricrogl_py = fullfile(tmpdir, 'mricrogl.py');
fid = fopen(mricrogl_py, 'wt');
fwrite(fid, mricrogl_command);
fclose(fid);

system(sprintf('/Applications/MRIcroGL.app/Contents/MacOS/MRIcroGL %s', mricrogl_py));

rmdir(tmpdir, 's');

end
