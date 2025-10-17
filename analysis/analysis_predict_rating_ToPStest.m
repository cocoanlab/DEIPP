function analysis_predict_rating_ToPStest(sj_num, topsdir, varargin)

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

addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));
addpath(genpath(fullfile(topsdir, 'functions', 'DCCcode')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
ratdir = fullfile(andir, 'rating', sj_id);
preddir = fullfile(andir, 'predict_rating', sj_id);

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

tasklabels = {'task-rating_run-1', 'task-rating_run-2', 'task-rating_run-3'};
func = [];
for task_i = 1:numel(tasklabels)
    funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabels{task_i}));
    funcmetafile = split(deblank(ls(funcmetafile)));
    for ses_i = 1:numel(funcmetafile)
        tempstruct = load(funcmetafile{ses_i}, 'func');
        func{task_i,ses_i} = change_bidsdir(tempstruct.func, bidsdir);
    end
end

TR = func{1,1}.raw.jsondat.RepetitionTime;
ndvol = ceil(30/TR);
yhdel = 13;

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'nbin', 'rbin_idx', 'tbin_idx');

topsparc_id = double(read_avw(fullfile(topsdir, 'data', 'atlas', 'Fan_et_al_atlas_r279_MNI_2mm.nii')));
topsparc_id = topsparc_id(:);
topsparc_id(topsparc_id == 0) = NaN;

topsname = fullfile(topsdir, 'model', 'ToPS_weight.mat');
if ~exist(topsname, 'file'); return; end
load(topsname, 'ToPS_w');

%% Apply ToPS

YToPS = cell(size(func));
[YToPS_rbin, YToPS_tbin] = deal(cell([size(func) numel(nbin)]));

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)

        fprintf('Working on run %d, ses %d...\n', task_i, ses_i);
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);
        
        func_dat = double(read_avw(func{task_i,ses_i}.den.bold));
        func_dat = reshape(func_dat, [], size(func_dat, 4));
        func_dat = func_dat(:, 1+ndvol+yhdel:end).';
        func_dat = func_dat(wh_keep, :);
        func_dat = splitapply(@(x) mean(x,2), func_dat, topsparc_id.');

        func_conn = DCC_jj(func_dat, 'whiten', 'simple');

        YToPS{task_i,ses_i} = sum(func_conn .* ToPS_w, 'omitnan').';

        for bin_i = 1:numel(nbin)

            YToPS_rbin{task_i,ses_i,bin_i} = splitapply(@(x) mean(x,1), YToPS{task_i, ses_i}, rbin_idx{task_i,ses_i,bin_i});
            YToPS_tbin{task_i,ses_i,bin_i} = splitapply(@(x) mean(x,1), YToPS{task_i, ses_i}, tbin_idx{task_i,ses_i,bin_i});

        end

    end
    
end

%% Save results

topstestname = fullfile(preddir, sprintf('%s_predict_rating_ToPS.mat', sj_id));
save(topstestname, 'YToPS', 'YToPS_rbin', 'YToPS_tbin');

fprintf('Done.\n');

end
