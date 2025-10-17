function analysis_predict_rating_trainsize(sj_num, varargin)

preptype = 'cocoan-preproc';
antype = 'cocoan-analysis';
bidsdir = fileparts(fileparts(mfilename('fullpath'))); % mfilename: bidsdir/code/~.m
if isempty(bidsdir); bidsdir = fileparts(pwd); end
parctype = 'indparc';

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
        end
    end
end

addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
parcdir = fullfile(andir, 'parcellation', sj_id);
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

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat_rbin', 'nbin', 'rbin_idx', 'tbin_idx');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s.mat', sj_id, parctype));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_rbin');

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

%% Prepare feature

func_conn = cell(size(func));

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)

        fprintf('Working on run %d, ses %d...\n', task_i, ses_i);
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);

        tempcf = ft_read_cifti_mod(func{task_i,ses_i}.cf.den.sbold);
        func_dat = tempcf.data;
        func_dat = func_dat(:, 1+ndvol+yhdel:end).';
        func_dat = func_dat(wh_keep, :);
        func_dat = splitapply(@(x) mean(x,2), func_dat, parc.id.');
        func_dat = func_dat(:, parc.whincl);
        
        func_conn{task_i,ses_i} = zscore(func_dat);
        func_conn{task_i,ses_i} = repmat(func_conn{task_i,ses_i}, 1, 1, size(func_conn{task_i,ses_i},2));
        func_conn{task_i,ses_i} = func_conn{task_i,ses_i} .* permute(func_conn{task_i,ses_i}, [1 3 2]);
        func_conn{task_i,ses_i} = func_conn{task_i,ses_i}(:, triu(true(size(func_conn{task_i,ses_i}, 2:3)), 1));
        
    end
    
end

%% Prepare data

X = func_conn_rbin(:,:,nbin==10);
Y = rating_dat_rbin(:,:,nbin==10);
idx_train = 1:size(Y,2);

rng('default');
trsz_nsub = 1:numel(idx_train)-2;
trsz_nrep = 100;
[trsz_allidx, trsz_cvidx] = deal(cell(1, numel(trsz_nsub)));
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

trsz_Yfit = cellfun(@(a) NaN(size(a,1), 1), func_conn, 'un', false);
trsz_Yfit = repmat(trsz_Yfit, 1, 1, trsz_nrep, numel(trsz_nsub));
[trsz_Yfit_rbin, trsz_Yfit_tbin] = deal(cell([size(trsz_Yfit) numel(nbin)]));

%% Training with changing data size

for sz_i = 1:numel(trsz_nsub)

    for rep_i = 1:trsz_nrep

        fprintf('Training size %d, Repitition %d ...\n', trsz_nsub(sz_i), rep_i);

        for k = idx_train

            wh_tr = setdiff(idx_train, k);
            wh_te = k;

            fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));

            wh_tr_sub = wh_tr(trsz_cvidx{sz_i}(:,rep_i,k));

            X_tr = cat(1, X{:,wh_tr_sub});
            Y_tr = cat(1, Y{:,wh_tr_sub});

            [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
            X_tr = X_tr*func_PC;
            [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', out.res.L_matlab(out.res.icv_bestidx(k)));
            regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
            bb = regfit_B;
            bb(regfit_B~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end)~=0)] \ Y_tr;
            
            w = [bb(1); func_PC * bb(2:end)];

            for task_i = 1:size(func,1)

                trsz_Yfit{task_i,k,rep_i,sz_i} = func_conn{task_i,k} * w(2:end) + w(1);

                for bin_i = 1:numel(nbin)

                    trsz_Yfit_rbin{task_i,k,rep_i,sz_i,bin_i} = splitapply(@(x) mean(x,1), ...
                        trsz_Yfit{task_i,k,rep_i,sz_i}, rbin_idx{task_i,k,bin_i});
                    trsz_Yfit_tbin{task_i,k,rep_i,sz_i,bin_i} = splitapply(@(x) mean(x,1), ...
                        trsz_Yfit{task_i,k,rep_i,sz_i}, tbin_idx{task_i,k,bin_i});

                end

            end

        end
    end

end

%% Calculate results and save model

trszname = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s.mat', sj_id, parctype));
save(trszname, 'trsz_nsub', 'trsz_nrep', 'trsz_cvidx', 'trsz_Yfit', 'trsz_Yfit_rbin', 'trsz_Yfit_tbin', '-v7.3');

fprintf('Done.\n');

end
