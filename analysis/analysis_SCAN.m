function analysis_SCAN(sj_num, mscgitdir, infomapdir, varargin)

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

setenv('DYLD_LIBRARY_PATH', '/opt/X11/lib/flat_namespace'); % for infomap
addpath(genpath(mscgitdir));
addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
mscdir = fullfile(andir, 'MSC', sj_id);
scandir = fullfile(andir, 'SCAN', sj_id);
if ~exist(scandir, 'dir'); mkdir(scandir); end

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

anatmetafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, 'ses-01', 'anat'));
load(anatmetafile, 'anat');
anat = change_bidsdir(anat, bidsdir);

tasklabel = 'task-rest';
func = [];
funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabel));
funcmetafile = split(deblank(ls(funcmetafile)));
for ses_i = 1:numel(funcmetafile)
    tempstruct = load(funcmetafile{ses_i}, 'func');
    func{ses_i} = change_bidsdir(tempstruct.func, bidsdir);
end

funcrefcf = ft_read_cifti_mod(func{1}.cf.den.sbold);
funcrefcf.data = [];
funcrefcf.pos(funcrefcf.brainstructure > 2, :) = [];
funcrefcf.brainstructure(funcrefcf.brainstructure > 2) = [];
funcrefcf.brainstructurelabel(3:end) = [];
funcrefcf = rmfield(funcrefcf,{'dim','transform'});
wh_nomw = funcrefcf.brainstructure ~= -1;

fsLRdir = fullfile(prepdir, 'ciftify', anat.basename, 'MNINonLinear', 'fsaverage_LR32k');

xdist = 30;
thresholds = [.0001 .0002 .0005 .001 .002 .005 .01 .02 .05];

%% Load data

% include only cortex
func_dat_all = [];
for ses_i = 1:numel(func)
    tempcf = ft_read_cifti_mod(func{ses_i}.cf.den.sbold);
    func_dat = double(tempcf.data(1:59412, logical(importdata(func{ses_i}.den.cen)))).';
    func_dat = func_dat - mean(func_dat);
    func_dat_all = [func_dat_all; func_dat];
end

func_corr = paircorr_mod(func_dat_all);
func_corr(isnan(func_corr)) = 0;
func_corr = FisherTransform(func_corr);

tmpdir = tempname(tempdir); mkdir(tmpdir);

surf_L = fullfile(fsLRdir, sprintf('%s.L.midthickness.32k_fs_LR.surf.gii', anat.basename));
surf_R = fullfile(fsLRdir, sprintf('%s.R.midthickness.32k_fs_LR.surf.gii', anat.basename));
system(sprintf('wb_command -surface-geodesic-distance-all-to-all %s %s', surf_L, fullfile(tmpdir, 'L.dconn.nii')));
system(sprintf('wb_command -surface-geodesic-distance-all-to-all %s %s', surf_R, fullfile(tmpdir, 'R.dconn.nii')));
distances = blkdiag(getfield(ft_read_cifti_mod(fullfile(tmpdir, 'L.dconn.nii')), 'data'), ...
    getfield(ft_read_cifti_mod(fullfile(tmpdir, 'R.dconn.nii')), 'data'));
distances(1:32492, 32493:64984) = 255;
distances(32493:64984, 1:32492) = 255;
distances = distances(wh_nomw, wh_nomw);

rmdir(tmpdir, 's');

%% Seed-based correlation along Lt. motor cortex

surf_L_gifti = gifti(surf_L);
BAfile = fullfile(fsLRdir, sprintf('%s.BA_exvivo.32k_fs_LR.dlabel.nii', anat.basename));
BAcf = ft_read_cifti_mod(BAfile);
LM_BA4p_idx = find(wh_nomw);
LM_BA4p_idx = LM_BA4p_idx(find(BAcf.data == 6));

switch sj_num
    case 1
        LM_seed_idx = [4461 3446 5225 19404 18755];
    case 2
        LM_seed_idx = [4461 3446 5225 19405 18752];
    case 3
        LM_seed_idx = [4461 3446 5279 19349 18752];
end

LM_edge_idx = reshape(cat(3, surf_L_gifti.faces, circshift(surf_L_gifti.faces, [0,1])), [], 2);
LM_edge_dist = sum((surf_L_gifti.vertices(LM_edge_idx(:,1), :).' - surf_L_gifti.vertices(LM_edge_idx(:,2), :).').^2).^0.5;
LM_graph = simplify(graph(LM_edge_idx(:,1), LM_edge_idx(:,2), LM_edge_dist));
LM_graph_BA4p = rmedge(LM_graph, find(~all(ismember(LM_graph.Edges.EndNodes, LM_BA4p_idx), 2)));
LM_seedpath_idx = NaN;
for p = 1:numel(LM_seed_idx)-1
    if all(ismember(LM_seed_idx(p:p+1), LM_BA4p_idx))
        LM_seedpath_idx = [LM_seedpath_idx(1:end-1) shortestpath(LM_graph_BA4p, LM_seed_idx(p), LM_seed_idx(p+1))];
    else
        LM_seedpath_idx = [LM_seedpath_idx(1:end-1) shortestpath(LM_graph, LM_seed_idx(p), LM_seed_idx(p+1))];
    end
end
[~, LM_seedpath_nomwidx] = ismember(LM_seedpath_idx, find(wh_nomw));

tempcf = funcrefcf;
tempcf.dimord = 'scalar_pos';
tempcf.data = func_corr(:, LM_seedpath_nomwidx);
tempcf.mapname = strcat({'Seed-based FC from vertex '}, cellstr(string(LM_seedpath_idx-1).'));
imgname = fullfile(scandir, sprintf('%s_SCAN_LM_seedcorr', sj_id));
ft_write_cifti_mod(imgname, tempcf);
wb_command = [];
for i = 1:size(tempcf.data, 2)
    wb_command = [wb_command, ...
        sprintf(['wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -column %d -pos-user %f %f -palette-name RBGYR20', ...
        ' -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_INSIDE %f %f -inversion POSITIVE_WITH_NEGATIVE;'], ...
        imgname, imgname, i, quantile(tempcf.data(:,i), [0.9 0.995 0.9 1]))];
end
system(wb_command);

%% Vertex-wise infomap with node-specific thresholding

tmpdir = tempname(tempdir); mkdir(tmpdir);

Run_Infomap_nodespecificthresh(func_corr, distances, xdist, thresholds, 0, tmpdir, [], 100, numel(thresholds), true, infomapdir);

communities = modify_clrfile('simplify',fullfile(tmpdir,'rawassn.txt'),10);
regularized = regularize(communities);
vertnetcf = funcrefcf;
vertnetcf.data = regularized;
ft_write_cifti_mod(fullfile(scandir,sprintf('%s_rawassn_minsize10_regularized',sj_id)), vertnetcf);

rmdir(tmpdir, 's');

%% Parcel-wise infomap with node-specific thresholding

parcels = ft_read_cifti_mod(fullfile(mscdir, sprintf('%s_parcels_edgethresh_0.5.dtseries.nii', sj_id)));
parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];

tcs = zeros(size(func_dat_all,1), length(parcelIDs));
for IDnum = 1:length(parcelIDs)
    tcs(:,IDnum) = mean(func_dat_all(:, parcels.data==parcelIDs(IDnum)),2);
end
parcel_corr = paircorr_mod(tcs);
parcel_corr(isnan(parcel_corr)) = 0;
parcel_corr = FisherTransform(parcel_corr);

parcel_centroids = zeros(length(parcelIDs),1);
for IDnum = 1:length(parcelIDs)
    parcelinds = find(parcels.data==parcelIDs(IDnum));
    within_parcel_distances = distances(parcelinds,parcelinds);
    [~,centroidind] = min(sum(within_parcel_distances,2));
    parcel_centroids(IDnum) = parcelinds(centroidind);
end
parcel_distances = distances(parcel_centroids,parcel_centroids);

tmpdir = tempname(tempdir); mkdir(tmpdir);

Run_Infomap_nodespecificthresh(parcel_corr, parcel_distances, xdist, thresholds, 0, tmpdir, [], 100, numel(thresholds), true, infomapdir);
communities = modify_clrfile('simplify',fullfile(tmpdir,'rawassn.txt'),1);
regularized = regularize(communities);

regularized_out = parcels;
regularized_out.data = zeros(size(parcels.data,1),size(regularized,2));
for IDnum = 1:length(parcelIDs)
    regularized_out.data(parcels.data==parcelIDs(IDnum),:) = repmat(regularized(IDnum,:),nnz(parcels.data==parcelIDs(IDnum)),1);
end
ft_write_cifti_mod(fullfile(scandir,sprintf('%s_rawassn_minsize1_regularized.dtseries.nii',sj_id)),regularized_out);

rmdir(tmpdir, 's');

%% Assigning SCAN labels

% Visual inspection and manual grouping of all inter-effector subnetworks
% (Gordon et al., 2023, Nature).
% For vertex-wise infomap, start with a 0.1% graph density
% (Gordon et al., 2020, PNAS). If the inter-effector subnetwork is not
% evident, proceed to sparser densities (0.05%, 0.02%, 0.01%).
% For parcel-wise Infomap, search across all thresholds to find
% subnetworks that match the vertex-wise infomap results.
%
% Vertex-wise infomap
% Participant 1
% 0.1%: 79 (lateral, middle), 50 (medial)
% Participant 2
% 0.1%: 39 (lateral, middle)
% 0.02%: 12 (medial)
% Participant 3
% 0.1%: 80 (lateral), 75(middle), 59 (medial)
%
% Parcel-wise infomap
% Participant 1
% 0.1%: 54, 56 (lateral, middle)
% 0.2%: 102 (medial)
% Participant 2
% 0.1%: 9 (lateral, middle), 21 (medial)
% Participant 2
% 0.2%: 5 (lateral, middle), 78 (medial)

tempcf = ft_read_cifti_mod(fullfile(scandir, sprintf('%s_rawassn_minsize10_regularized.dtseries.nii', sj_id)));
regularized = tempcf.data;
tempcf = ft_read_cifti_mod(fullfile(mscdir, sprintf('%s_rawassn_minsize400_regularized_recolored.dscalar.nii', sj_id)));
switch sj_num
    case 1
        tempcf.data(ismember(regularized(:, thresholds==.001), [50 79])) = 1.5;
    case 2
        tempcf.data(ismember(regularized(:, thresholds==.001), 39)) = 1.5;
        tempcf.data(ismember(regularized(:, thresholds==.0002), 12)) = 1.5;
    case 3
        tempcf.data(ismember(regularized(:, thresholds==.001), [59 75 80])) = 1.5;
end
imgname = fullfile(scandir, sprintf('%s_rawassn_minsize400_regularized_recolored_addSCAN', sj_id));
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));

tempcf = ft_read_cifti_mod(fullfile(scandir, sprintf('%s_rawassn_minsize1_regularized.dtseries.nii', sj_id)));
regularized = tempcf.data;
tempcf = ft_read_cifti_mod(fullfile(mscdir, sprintf('%s_rawassn_minsize4_regularized_recolored.dscalar.nii', sj_id)));
switch sj_num
    case 1
        tempcf.data(ismember(regularized(:, thresholds==.001), [54 56])) = 1.5;
        tempcf.data(ismember(regularized(:, thresholds==.002), 102)) = 1.5;
    case 2
        tempcf.data(ismember(regularized(:, thresholds==.001), [9 21])) = 1.5;
    case 3
        tempcf.data(ismember(regularized(:, thresholds==.002), [5 78])) = 1.5;
end
parcelnets = zeros(length(parcelIDs),1);
for IDnum = 1:length(parcelIDs)
    parcelnets(IDnum) = unique(tempcf.data(parcels.data==parcelIDs(IDnum)));
end
imgname = fullfile(scandir, sprintf('%s_rawassn_minsize4_regularized_recolored_addSCAN', sj_id));
ft_write_cifti_mod(imgname, tempcf);
system(sprintf('wb_command -cifti-palette %s.dscalar.nii MODE_USER_SCALE %s.dscalar.nii -pos-user 1 18 -palette-name power_surf', imgname, imgname));
dlmwrite([imgname '.txt'], parcelnets, 'delimiter', ' ');

end
