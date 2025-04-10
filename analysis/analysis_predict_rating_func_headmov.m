function analysis_predict_rating_func_headmov(sj_num, varargin)

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

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
preddir = fullfile(andir, 'predict_rating_func', sj_id);

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

tasklabel = {'task-rating_run-1', 'task-rating_run-2', 'task-rating_run-3'};
func = [];
for task_i = 1:numel(tasklabel)
    funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabel{task_i}));
    funcmetafile = split(deblank(ls(funcmetafile)));
    for ses_i = 1:numel(funcmetafile)
        tempstruct = load(funcmetafile{ses_i}, 'func');
        func{task_i,ses_i} = change_bidsdir(tempstruct.func, bidsdir);
    end
end

TR = func{1,1}.raw.jsondat.RepetitionTime;
nvol = func{1,1}.rmd.nvol;
ndvol = ceil(30/TR);
yhdel = 13;

%% Load data

obsname = fullfile(preddir, sprintf('%s_predict_rating_observation.mat', sj_id));
if ~exist(obsname, 'file'); return; end
load(obsname, 'nbin', 'rbin_idx', 'tbin_idx');

%% Head motion data

fd_nocens = {};
fd = {};

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);
        
        rp = importdata(func{task_i,ses_i}.mc.rp);
        rp = [rp(:,4:6), rp(:,1:3) * pi / 180 * 50];
        fd{task_i,ses_i} = sum(abs(diff(rp([1 1:end], :))), 2);
        fd{task_i,ses_i} = fd{task_i,ses_i}(1+ndvol+yhdel:end);
        fd_nocens{task_i,ses_i} = fd{task_i,ses_i};
        fd{task_i,ses_i} = fd{task_i,ses_i}(wh_keep);

    end
    
end

[fd_rbin, fd_tbin] = deal(cell([size(func) numel(nbin)]));

for bin_i = 1:numel(nbin)
    fd_rbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), fd, rbin_idx(:,:,bin_i), 'un', false);
    fd_tbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), fd, tbin_idx(:,:,bin_i), 'un', false);
end

%% Save results

headmovname = fullfile(preddir, sprintf('%s_predict_rating_headmov.mat', sj_id));
save(headmovname, 'fd_nocens', 'fd', 'fd_rbin', 'fd_tbin');

fprintf('Done.\n');

end
