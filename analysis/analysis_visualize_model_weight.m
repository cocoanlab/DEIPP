function analysis_visualize_model_weight(sj_num, spmdir, mricrogldir, varargin)

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

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

figbasename = fullfile(preddir, sprintf('%s_predict_rating_weight_%s', sj_id, parctype));

%% Prepare model weight table

w = refmt_r(out.beta_best_edge.all(2:end));
wpos = sum(w .* double(w > 0), 2);
wneg = sum(w .* double(w < 0), 2);

w_tbl = table((1:numel(parc.net)).', string(parc.names), round(parc.medpos), string(parc.netnames(parc.net)), ...
    'VariableNames', {'ROI', 'name', 'MNI', 'netnames'});
w_tbl = w_tbl(parc.whincl, :);
w_tbl = addvars(w_tbl, wpos, wneg, 'After', 1, 'NewVariableNames', {'wpos', 'wneg'});
head(sortrows(w_tbl, 'wpos', 'descend'), 20);
head(sortrows(w_tbl, 'wneg', 'ascend'), 20);

switch sj_num
    case 1
        wpos_cl = [0 0.07];
        wneg_cl = [-0.06 0];
    case 2
        wpos_cl = [0 0.11];
        wneg_cl = [-0.14 0];
end

%% Visualize cortex

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = 1:numel(w_tbl.ROI)
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wpos(i);
end
tempcf.mapname = {'Sum of positive weights'};
imgname = sprintf('%s_pos', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user %f %f -palette-name red-yellow', imgname, imgname, wpos_cl));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = 1:numel(w_tbl.ROI)
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wneg(i);
end
tempcf.mapname = {'Sum of negative weights'};
imgname = sprintf('%s_neg', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -neg-user %f %f -palette-name blue-lightblue -inversion POSITIVE_WITH_NEGATIVE', imgname, imgname, flip(wneg_cl)));

%% Visualize cerebellum

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
cbvert = NaN(size(cbparc.cdata));
cb_idx = sort(w_tbl.ROI(strcmp(w_tbl.netnames, 'CB')), 'ascend');
for i = 1:numel(cb_idx)
    cbvert(cbparc.cdata == i) = w_tbl.wpos(w_tbl.ROI == cb_idx(i));
end
cbvert(cbvert == 0) = NaN;
wposcol = cols.fx_interp(cols.wb_ry, 1000);
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', wposcol, 'cscale', wpos_cl, 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_pos_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
cbvert = NaN(size(cbparc.cdata));
cb_idx = sort(w_tbl.ROI(strcmp(w_tbl.netnames, 'CB')), 'ascend');
for i = 1:numel(cb_idx)
    cbvert(cbparc.cdata == i) = w_tbl.wneg(w_tbl.ROI == cb_idx(i));
end
cbvert(cbvert == 0) = NaN;
wnegcol = cols.fx_interp(flipud(cols.wb_blb), 1000);
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', wnegcol, 'cscale', wneg_cl, 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_neg_cerebellum.png', figbasename);
exportgraphics(gcf, figname);
close all;

%% Visualize subcortex

tmpdir = tempname(tempdir); mkdir(tmpdir);

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(w_tbl.netnames, {'THA', 'BG'})).'
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wpos(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'wpos_thabg');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(w_tbl.netnames, {'HCAMY', 'BS'})).'
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wpos(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'wpos_hcamybs');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(w_tbl.netnames, {'THA', 'BG'})).'
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wneg(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'wneg_thabg');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
for i = find(ismember(w_tbl.netnames, {'HCAMY', 'BS'})).'
    tempcf.data(parc.id == w_tbl.ROI(i)) = w_tbl.wneg(i);
end
tempcf.mapname = {''};
imgname = fullfile(tmpdir, 'wneg_hcamybs');
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-separate %s.dscalar.nii COLUMN -volume-all %s.nii.gz', imgname, imgname));

make_CLUT(cols.wb_ry, fullfile(mricroglclutdir, 'wb_ry.clut'));
make_CLUT(flipud(cols.wb_blb), fullfile(mricroglclutdir, 'wb_blb.clut'));

mricrogl_command = ['import gl', newline, ...
    'gl.resetdefaults()', newline, ...
    'gl.backcolor(255,255,255)', newline, ...
    'gl.colorbarposition(0)', newline, ...
    sprintf('gl.loadimage("%s/data/standard/MNI152_T1_1mm_brain.nii.gz")', getenv('FSLDIR')), newline, ...
    'gl.overlayloadsmooth(0)', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'wpos_thabg')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', wpos_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_ry")', newline, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_pos_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'wpos_hcamybs')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', wpos_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_ry")', newline, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_pos_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'wneg_thabg')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', wneg_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_blb")', newline, ...
    'gl.mosaic("A 18 9 -1 H 0.42")', newline, ...
    sprintf('gl.savebmp("%s_neg_z_18_9_-1.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, ...
    sprintf('gl.overlayload("%s.nii.gz")', fullfile(tmpdir, 'wneg_hcamybs')), newline, ...
    sprintf('gl.minmax(1, %.4f, %.4f)', wneg_cl), newline, ...
    'gl.zerointensityinvisible(1, 1)', newline, ...
    'gl.colorname(1, "wb_blb")', newline, ...
    'gl.mosaic("Z -25 -2 25 H 0.48")', newline, ...
    sprintf('gl.savebmp("%s_neg_x_-25_-2_25.png")', figbasename), newline, ...
    'gl.overlaycloseall()', newline];

mricrogl_command = [mricrogl_command, 'gl.quit()'];

mricrogl_py = fullfile(tmpdir, 'mricrogl.py');
fid = fopen(mricrogl_py, 'wt');
fwrite(fid, mricrogl_command);
fclose(fid);

system(sprintf('%s %s', mricrogldir, mricrogl_py));

rmdir(tmpdir, 's');

delete(fullfile(mricroglclutdir, 'wb_ry.clut'));
delete(fullfile(mricroglclutdir, 'wb_blb.clut'));

end
