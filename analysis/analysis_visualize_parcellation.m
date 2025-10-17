function analysis_visualize_parcellation(sj_num, spmdir, mricrogldir, varargin)

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

cols = call_colors;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

figbasename = fullfile(parcdir, sprintf('%s_parcellation_%s', sj_id, parctype));

%% Visualize cortex

tempcf = parc.cf;
tempcf.data = parc.id;
tempcf.mapname = {'ROI number'};
imgname = sprintf('%s_roi', figbasename);
ft_write_cifti_mod(imgname, tempcf);

tempcf = parc.cf;
tempcf.data = parc.id .* 0 + 1;
tempcf.mapname = {'ROI binary'};
imgname = sprintf('%s_binaryroi', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 2 -palette-name Gray_Interp_Positive', imgname, imgname));

tempcf = parc.cf;
tempcf.data = NaN(size(parc.id));
tempcf.data(~isnan(parc.id)) = parc.netorig(parc.id(~isnan(parc.id)));
tempcf.mapname = {'ROI network label'};
imgname = sprintf('%s_netroi', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));

tempcf = parc.cf;
tempcf.data = double(ismember(parc.id, find(~parc.whincl)));
tempcf.mapname = {'ROI in the visual network and occipital areas, excluded'};
imgname = sprintf('%s_visroi', figbasename);
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 0 10000 -palette-name spectral', imgname, imgname));

%% Visualize cerebellum

cbparc = gifti(fullfile(spmdir, 'toolbox/suit/flatmap/MDTB_10Regions.label.gii'));
cbparc.cdata = double(cbparc.cdata);
cbvert = NaN(size(cbparc.cdata));
roi_idx = 1:max(cbparc.cdata);
for j = 1:numel(roi_idx)
    cbvert(cbparc.cdata == roi_idx(j)) = 1 - (j-1)/numel(roi_idx);
end
cbvert(cbvert == 0) = NaN;
netcol = lab2rgb(rgb2lab(cols.wb_power(ismember(parc.netnames, 'CB'), :)) + [linspace(-3, 3, 255).' zeros(255, 2)]);
suit_plotflatmap(cbvert, 'type', 'func', 'cmap', netcol, 'cscale', [0 1], 'border', []);
axis off; set(gcf, 'Color', 'w');
figname = sprintf('%s_netroi_cerebellum.png', figbasename);
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
    sprintf('gl.savebmp("%s_z_18_9_-1.png")', sprintf('%s_netroi', figbasename)), newline, ...
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
    sprintf('gl.savebmp("%s_x_-25_-2_25.png")', sprintf('%s_netroi', figbasename)), newline, ...
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
