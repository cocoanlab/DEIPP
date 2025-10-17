function analysis_MSC(sj_num, mscgitdir, infomapdir, varargin)

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
if ~exist(mscdir, 'dir'); mkdir(mscdir); end

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

edgethresh = .5;
xdist = 30;
thresholds = [.003 .004 .005:.005:.05];

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

%% Vertex-wise infomap

tmpdir = tempname(tempdir); mkdir(tmpdir);

Run_Infomap_2015(func_corr, distances, xdist, thresholds, 0, tmpdir, numel(thresholds), [], infomapdir);

communities = modify_clrfile('simplify',fullfile(tmpdir,'rawassn.txt'),400);
regularized = regularize(communities);
vertnetcf = funcrefcf;
vertnetcf.data = regularized;
ft_write_cifti_mod(fullfile(mscdir,sprintf('%s_rawassn_minsize400_regularized',sj_id)), vertnetcf);
% min 30 mm^2 option was added following the Gordon et al., 2017
consensus_maker_knowncolors(fullfile(mscdir,sprintf('%s_rawassn_minsize400_regularized.dtseries.nii',sj_id)), [], which('Networks_template_cleaned.dscalar.nii'), 0, [], 30);
make_block_diagram(fullfile(mscdir,sprintf('%s_rawassn_minsize400_regularized_allcolumns_recolored.dscalar.nii',sj_id)), thresholds);
close all;
cifti_to_border_v2(fullfile(mscdir,sprintf('%s_rawassn_minsize400_regularized_recolored.dscalar.nii',sj_id)), 1, 1, 'default');

rmdir(tmpdir, 's');

%% Individualized parcellation

funccorrcf = funcrefcf;
funccorrcf.data = func_corr;

tmpdir = tempname(tempdir); mkdir(tmpdir);

currdir = pwd;
surface_parcellation_singlesub(anat.basename, funccorrcf, fsLRdir, 1, 0, tmpdir, round(maxNumCompThreads*0.8));
parcel_creator_cifti(fullfile(tmpdir, 'corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii'), ...
    fullfile(tmpdir, sprintf('%s_parcels', sj_id)), edgethresh);
cd(currdir);
funccorrcf = [];

copyfile(fullfile(tmpdir, '*'), mscdir);
rmdir(tmpdir, 's');

%% Parcel-wise infomap

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

Run_Infomap_2015(parcel_corr, parcel_distances, xdist, thresholds, 0, tmpdir, numel(thresholds), [], infomapdir);
communities = modify_clrfile('simplify',fullfile(tmpdir,'rawassn.txt'),4);
regularized = regularize(communities);

regularized_out = parcels;
regularized_out.data = zeros(size(parcels.data,1),size(regularized,2));
for IDnum = 1:length(parcelIDs)
    regularized_out.data(parcels.data==parcelIDs(IDnum),:) = repmat(regularized(IDnum,:),nnz(parcels.data==parcelIDs(IDnum)),1);
end
ft_write_cifti_mod(fullfile(mscdir,sprintf('%s_rawassn_minsize4_regularized.dtseries.nii',sj_id)),regularized_out);
consensus_maker_knowncolors(fullfile(mscdir,sprintf('%s_rawassn_minsize4_regularized.dtseries.nii',sj_id)),[],...
    fullfile(mscdir,sprintf('%s_rawassn_minsize400_regularized_recolored.dscalar.nii',sj_id)),0,[],[],fullfile(mscdir, sprintf('%s_parcels_edgethresh_0.5.dtseries.nii', sj_id)))

rmdir(tmpdir, 's');

%% Spring network

spring_embedding_func_easy_crossthresh(parcel_corr,fullfile(mscdir,sprintf('%s_rawassn_minsize4_regularized_recolored.txt',sj_id)),...
    1,25,parcel_distances,xdist,thresholds,fullfile(mscdir,sprintf('%s_spring_embed',sj_id)));
close all;

end
