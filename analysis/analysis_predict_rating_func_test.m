function analysis_predict_rating_func_test(sj_num, opt, varargin)

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

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

obsname = fullfile(preddir, sprintf('%s_predict_rating_observation.mat', sj_id));
if ~exist(obsname, 'file'); return; end
load(obsname, 'nbin', 'rbin_idx', 'tbin_idx');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s_%s.mat', ...
    sj_id, opt.dfc, opt.parc));
if ~exist(featname, 'file'); return; end
load(featname, 'roi_id', 'wh_exclude');

%% Prepare feature

func_dat = {};

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);
        
        tempcf = ft_read_cifti_mod(func{task_i,ses_i}.cf.bold);
        func_dat{task_i,ses_i} = tempcf.data;
        func_dat{task_i,ses_i} = func_dat{task_i,ses_i}(:, 1+ndvol+yhdel:end).';
        func_dat{task_i,ses_i} = func_dat{task_i,ses_i}(wh_keep, :);
        func_dat{task_i,ses_i} = splitapply(@(x) mean(x,2), func_dat{task_i,ses_i}, roi_id.');
        func_dat{task_i,ses_i} = func_dat{task_i,ses_i}(:,wh_exclude);
        
    end
    
end

switch opt.dfc
    case 'ets'
        func_conn = cellfun(@(a) zscore(a), func_dat, 'un', false);
        func_conn = cellfun(@(a) repmat(a,1,1,size(a,2)), func_conn, 'un', false);
        func_conn = cellfun(@(a) a .* permute(a, [1 3 2]), func_conn, 'un', false);
        func_conn = cellfun(@(a) a(:,triu(true(size(a,2:3)),1)), func_conn, 'un', false);
    case {'sw10', 'sw20', 'sw30'}
        switch opt.dfc
            case 'sw10'
                nwin = round(10/TR);
            case 'sw20'
                nwin = round(20/TR);
            case 'sw30'
                nwin = round(30/TR);
        end
        func_conn = cellfun(@(a) NaN(size(a,1), nchoosek(size(a,2),2)), func_dat, 'un', false);
        for i = 1:numel(func_dat)
            for j = 1:size(func_dat{i},1)
                targ_idx = max(j-floor(nwin/2), 1) : min(j+ceil(nwin/2)-1, size(func_dat{i},1));
                func_conn{i}(j,:) = atanh(refmt_r(corr(func_dat{i}(targ_idx, :))));
            end
        end
end

%% Apply model

Yfit = cell(size(func));
[Yfit_rbin, Yfit_tbin] = deal(cell([size(func) numel(nbin)]));

for j = 1:size(func_conn, 2)
    if ismember(j, out.idx_hout)
        w = out.beta_best_edge.all;
    else
        w = out.beta_best_edge.ocv{j};
    end
    Yfit(:,j) = cellfun(@(a) a * w(2:end) + w(1), func_conn(:,j), 'un', false);
end

for bin_i = 1:numel(nbin)
    Yfit_rbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), Yfit, rbin_idx(:,:,bin_i), 'un', false);
    Yfit_tbin(:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), Yfit, tbin_idx(:,:,bin_i), 'un', false);
end

%% Save results

testname = fullfile(preddir, sprintf('%s_predict_rating_test_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
save(testname, 'Yfit', 'Yfit_rbin', 'Yfit_tbin');

fprintf('Done.\n');

end