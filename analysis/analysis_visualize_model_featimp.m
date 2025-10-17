function analysis_visualize_model_featimp(sj_num, spmdir, mricrogldir, varargin)

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

figbasename = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s', sj_id, parctype));

%% Prepare feature importance table

fimp_tbl = table((1:numel(parc.net)).', string(parc.names), round(parc.medpos), string(parc.netnames(parc.net)), ...
    'VariableNames', {'ROI', 'name', 'MNI', 'netnames'});
fimp_tbl = fimp_tbl(parc.whincl, :);
fimp_tbl = addvars(fimp_tbl, fimp.reg, 'After', 1, 'NewVariableNames', 'imp');
head(sortrows(fimp_tbl, 'imp', 'descend'), 20);

switch sj_num
    case 1
        fimp_cl = [0 0.05];
    case 2
        fimp_cl = [0 0.02];
end

%% Top five, dot

topk = 5;
fimp_topk = maxk(fimp_tbl.imp, topk);
dotcols = interp1(linspace(fimp_cl(1), fimp_cl(2), size(cols.wb_cw_w,1)), cols.wb_cw_w, min(fimp_topk, fimp_cl(2)));
hold on;
for i = 1:topk
    line([0 fimp_topk(i)], [6-i 6-i], 'Color', dotcols(i,:), 'LineWidth', 1.5);
    scatter(fimp_topk(i), 6-i, 100, dotcols(i,:), 'filled');
end
set(gca, 'LineWidth', 2, 'YLim', [0.65 topk+0.35], 'YTick', [], ...
    'TickLength', [0.015 0.015], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);
switch sj_num
    case 1
        set(gca, 'XLim', [0 0.053], 'XTick', 0:0.02:0.04, 'XGrid', 'on', 'GridLineStyle', '--');
    case 2
        set(gca, 'XLim', [0 0.025], 'XTick', 0:0.01:0.02, 'XGrid', 'on', 'GridLineStyle', '--');
end

resize_axes(gca, [130 398]);
set(gcf, 'color', 'w');

figname = sprintf('%s_topfive_dot.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname);
close all;

%% Visualize cortex

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = 1:numel(fimp_tbl.ROI)
    tempcf.data(parc.id == fimp_tbl.ROI(i)) = fimp_tbl.imp(i);
end
tempcf.mapname = {'Importance'};
imgname = sprintf('%s', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user %f %f -neg-user %f %f -palette-name cool-warm', imgname, imgname, [fimp_cl -1 -1]));

%% Visualize cerebellum

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
cbvert = NaN(size(cbparc.cdata));
cb_idx = sort(fimp_tbl.ROI(strcmp(fimp_tbl.netnames, 'CB')), 'ascend');
for i = 1:numel(cb_idx)
    cbvert(cbparc.cdata == i) = fimp_tbl.imp(fimp_tbl.ROI == cb_idx(i));
end
cbvert(cbvert == 0) = NaN;
fimpcol = cols.fx_interp(cols.wb_cw_w, 1000);
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', fimpcol, 'cscale', fimp_cl, 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

%% Visualize subcortex

tmpdir = tempname(tempdir); mkdir(tmpdir);

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(fimp_tbl.netnames, {'THA', 'BG'})).'
    tempcf.data(parc.id == fimp_tbl.ROI(i)) = fimp_tbl.imp(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'fimp_thabg');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(fimp_tbl.netnames, {'HCAMY', 'BS'})).'
    tempcf.data(parc.id == fimp_tbl.ROI(i)) = fimp_tbl.imp(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'fimp_hcamybs');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

make_CLUT(cols.wb_cw_w, fullfile(mricroglclutdir, 'wb_cw_w.clut'));

mricrogl_command = ['import gl', newline, ...
    'gl.resetdefaults()', newline, ...
    'gl.backcolor(255,255,255)', newline, ...
    'gl.colorbarposition(0)', newline, ...
    sprintf('gl.loadimage("%s/data/standard/MNI152_T1_1mm_brain.nii.gz")', getenv('FSLDIR')), newline, ...
    'gl.overlayloadsmooth(0)', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'fimp_thabg')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', fimp_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_cw_w")', newline, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'fimp_hcamybs')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', fimp_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_cw_w")', newline, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, 'gl.quit()'];

mricrogl_py = fullfile(tmpdir, 'mricrogl.py');
fid = fopen(mricrogl_py, 'wt');
fwrite(fid, mricrogl_command);
fclose(fid);

system(sprintf('%s %s', mricrogldir, mricrogl_py));

rmdir(tmpdir, 's');

delete(fullfile(mricroglclutdir, 'wb_cw_w.clut'));

fprintf('Done.\n');

end
