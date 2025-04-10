function analysis_predict_rating_func_observation(sj_num, varargin)

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

nbin = [2 5 10];

%% Prepare observation

rating_dat_nocens = {};
rating_dat = {};
rbin_idx = {};
rbin_thr_bad = 5;
tbin_idx = {};

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)
        
        ses_id = sprintf('ses-%.2d', ses_i);
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);

        temptsv = fullfile(bidsdir, sj_id, ses_id, 'func', sprintf('%s_%s_%s_resp.tsv', sj_id, ses_id, tasklabel{task_i}));
        rating_dat{task_i,ses_i} = importdata(temptsv);
        rating_dat{task_i,ses_i} = interp1(rating_dat{task_i,ses_i}(:,1), ...
            rating_dat{task_i,ses_i}(:,2), (22*TR-10+0.5*TR : TR : 22*TR-10+0.5*TR+(nvol-1)*TR).');
        rating_dat{task_i,ses_i} = rating_dat{task_i,ses_i}(1+ndvol:end-yhdel);
        rating_dat_nocens{task_i,ses_i} = rating_dat{task_i,ses_i};
        rating_dat{task_i,ses_i} = rating_dat{task_i,ses_i}(wh_keep);
        
        for bin_i = 1:numel(nbin)
            
            rating_dscrank = rating_dat{task_i,ses_i};
            [~, wh_sortbin] = sort(rating_dscrank, 'ascend');
            rating_dscrank(wh_sortbin) = 1:numel(rating_dscrank);
            rating_quantile = sum(rating_dscrank > [-Inf quantile(rating_dscrank, (1:nbin(bin_i)-1)./nbin(bin_i))], 2);
%             rating_quantile = rating_dat{task_i,ses_i};
%             rating_quantile = sum(rating_quantile > [-Inf quantile(rating_quantile, nbin(bin_i)-1)], 2);
            wh_quantile_diff = [1; find(diff(rating_quantile([1 1:end])))];
            wh_quantile_diff_bad = find(diff(wh_quantile_diff) < rbin_thr_bad);
            wh_quantile_diff_bad = [wh_quantile_diff(wh_quantile_diff_bad) wh_quantile_diff(wh_quantile_diff_bad+1)];
            wh_quantile_bad = false(numel(rating_quantile),1);
            for diff_i = 1:size(wh_quantile_diff_bad,1)
                wh_quantile_bad(wh_quantile_diff_bad(diff_i,1) : wh_quantile_diff_bad(diff_i,2)-1) = true;
            end
            rbin_idx{task_i,ses_i,bin_i} = rating_quantile;
            rbin_idx{task_i,ses_i,bin_i}(wh_quantile_bad) = NaN;

            tbin_idx{task_i,ses_i,bin_i} = discretize(1:numel(rating_dat{task_i,ses_i}), nbin(bin_i)).';
            
        end
        
    end
    
end

[rating_dat_rbin, rating_dat_tbin] = deal(cell([size(func) numel(nbin)]));

for bin_i = 1:numel(nbin)
    rating_dat_rbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), rating_dat, rbin_idx(:,:,bin_i), 'un', false);
    rating_dat_tbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), rating_dat, tbin_idx(:,:,bin_i), 'un', false);
end

%% Save observation

obsname = fullfile(preddir, sprintf('%s_predict_rating_observation.mat', sj_id));
save(obsname, 'rating_dat_nocens', 'rating_dat', 'rating_dat_rbin', 'rating_dat_tbin', ...
    'nbin', 'rbin_idx', 'tbin_idx');

fprintf('Done.\n');

end

