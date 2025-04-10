function analysis_predict_rating_func_trainsize(sj_num, opt, varargin)

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
load(obsname, 'rating_dat_rbin', 'rating_dat_tbin', 'nbin', 'rbin_idx', 'tbin_idx');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s_%s.mat', ...
    sj_id, opt.dfc, opt.parc));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_rbin', 'func_conn_tbin', 'roi_id', 'wh_exclude');

%% Prepare feature

func_dat = {};

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)
        
        ses_id = sprintf('ses-%.2d', ses_i);
        
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

%% Prepare data

switch opt.bintype
    case 'rbin'
        X = func_conn_rbin(:,:,nbin==opt.nbin);
        Y = rating_dat_rbin(:,:,nbin==opt.nbin);
    case 'tbin'
        X = func_conn_tbin(:,:,nbin==opt.nbin);
        Y = rating_dat_tbin(:,:,nbin==opt.nbin);
end

idx_hout = out.idx_hout;
idx_train = out.idx_train;

rng('default');
trsz_nsub = 1:numel(idx_train)-2;
trsz_nrep = 100;
trsz_allidx = [];
trsz_cvidx = [];
for sz_i = 1:numel(trsz_nsub)
    trsz_allidx{sz_i} = NaN(trsz_nsub(sz_i), trsz_nrep);
    trsz_cvidx{sz_i} = NaN(trsz_nsub(sz_i), trsz_nrep, numel(idx_train));
    for rep_i = 1:trsz_nrep
        trsz_allidx{sz_i}(:,rep_i) = randperm(numel(idx_train), trsz_nsub(sz_i));
        for k = 1:numel(idx_train)
            trsz_cvidx{sz_i}(:,rep_i,k) = randperm(numel(idx_train)-1, trsz_nsub(sz_i));
        end
    end
end

trsz_Yfit = cellfun(@(a) NaN(size(a,1), 1), func_dat, 'un', false);
trsz_Yfit = repmat(trsz_Yfit, 1, 1, trsz_nrep, numel(trsz_nsub));
[trsz_Yfit_rbin, trsz_Yfit_tbin] = deal(cell([size(trsz_Yfit) numel(nbin)]));

%% Training with changing data size

for sz_i = 1:numel(trsz_nsub)

    for rep_i = 1:trsz_nrep

        fprintf('Training size %d, Repitition %d ...\n', trsz_nsub(sz_i), rep_i);

        [beta, beta_edge] = deal([]);
        [beta.all, beta_edge.all] = deal(cell(1, 1));
        [beta.ocv, beta_edge.ocv] = deal(cell(size(Y,2), 1));

        % Training: All

        wh_tr = idx_train;
        wh_te = idx_hout;

        if ~isempty(wh_te)

            fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));

            wh_tr_sub = wh_tr(trsz_allidx{sz_i}(:,rep_i));

            X_tr = cat(1, X{:,wh_tr_sub});
            Y_tr = cat(1, Y{:,wh_tr_sub});

            switch opt.algorithm
                case {'LASSOPCR', 'LASSOOLSPCR'}
                    [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
                    X_tr = X_tr*func_PC;
                    [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', out.res.L_matlab(out.res.ocv_bestidx));
                    regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
                    beta.all = regfit_B;
                    if strcmp(opt.algorithm, 'LASSOOLSPCR')
                        beta.all(regfit_B~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end)~=0)] \ Y_tr;
                    end
                    beta_edge.all = [beta.all(1); func_PC * beta.all(2:end)];
            end

        end

        % Training: CV

        for k = idx_train

            wh_tr = setdiff(idx_train, k);
            wh_te = k;

            fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));

            wh_tr_sub = wh_tr(trsz_cvidx{sz_i}(:,rep_i,k));

            X_tr = cat(1, X{:,wh_tr_sub});
            Y_tr = cat(1, Y{:,wh_tr_sub});

            switch opt.algorithm
                case {'LASSOPCR', 'LASSOOLSPCR'}
                    [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
                    X_tr = X_tr*func_PC;
                    [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', out.res.L_matlab(out.res.icv_bestidx(k)));
                    regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
                    beta.ocv{k} = regfit_B;
                    if strcmp(opt.algorithm, 'LASSOOLSPCR')
                        beta.ocv{k}(regfit_B~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end)~=0)] \ Y_tr;
                    end
                    beta_edge.ocv{k} = [beta.ocv{k}(1); func_PC * beta.ocv{k}(2:end)];
            end

        end
        
        for j = 1:size(func_conn, 2)
            if ismember(j, out.idx_hout)
                w = beta_edge.all;
            else
                w = beta_edge.ocv{j};
            end
            trsz_Yfit(:,j,rep_i,sz_i) = cellfun(@(a) a * w(2:end) + w(1), func_conn(:,j), 'un', false);
        end
    end

end

for bin_i = 1:numel(nbin)
    trsz_Yfit_rbin(:,:,:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), trsz_Yfit, ...
        repmat(rbin_idx(:,:,bin_i), 1, 1, trsz_nrep, numel(trsz_nsub)), 'un', false);
    trsz_Yfit_tbin(:,:,:,:,bin_i) = cellfun(@(a,b) ...
        splitapply(@(aa) mean(aa,1), a, b), trsz_Yfit, ...
        repmat(tbin_idx(:,:,bin_i), 1, 1, trsz_nrep, numel(trsz_nsub)), 'un', false);
end

%% Calculate results and save model

trszname = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
save(trszname, 'trsz_nsub', 'trsz_nrep', 'trsz_allidx', 'trsz_cvidx', 'trsz_Yfit', 'trsz_Yfit_rbin', 'trsz_Yfit_tbin', '-v7.3');

fprintf('Done.\n');

end