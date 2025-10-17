function analysis_predict_rating_featimp(sj_num, varargin)

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

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat_tbin', 'nbin');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s.mat', sj_id, parctype));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_tbin');

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

%% Prepare feature importance

X = func_conn_tbin(:,:,nbin==10);
Y = rating_dat_tbin(:,:,nbin==10);

w = out.beta_best_edge.ocv.';

Y_cat = cat(1, Y{:});

w_int = cellfun(@(a) a(1, :), w);
w_beta = cellfun(@(a) a(2:end, :), w, 'un', false);

X_dot_w = cellfun(@(a,b) a .* b.', X, repmat(w_beta, size(X,1), 1), 'un', false);
X_dot_w = cat(1, X_dot_w{:});
X_dotprod_w_orig = sum(X_dot_w, 2) + reshape(repmat(w_int, size(X{1,1},1)*size(X,1), 1), [], 1);
corr_orig = corr(X_dotprod_w_orig, Y_cat);

rng('default');
nrep = 10000;
perm_idx = NaN(numel(Y_cat), nrep);
for rep_i = 1:nrep
    perm_idx_each = reshape(1:size(perm_idx,1), [], size(Y,1), size(Y,2));
    perm_idx_each = perm_idx_each(:, :, randperm(size(perm_idx_each,3)));
    for ses_i = 1:size(perm_idx_each,3)
        perm_idx_each(:,:,ses_i) = perm_idx_each(:, randperm(size(perm_idx_each,2)), ses_i);
        for task_i = 1:size(perm_idx_each,2)
            perm_idx_each(:,task_i,ses_i) = perm_idx_each(randperm(size(perm_idx_each,1)), task_i, ses_i);
        end
    end
    perm_idx(:,rep_i) = perm_idx_each(:);
end

%% Region-wise feature importance

X_dot_w_reg = cellfun(@refmt_r, num2cell(X_dot_w, 2), 'un', false);
X_dot_w_reg = cellfun(@sum, X_dot_w_reg, 'un', false);
X_dot_w_reg = cat(1, X_dot_w_reg{:});
X_dotprod_w_reg_rem = X_dotprod_w_orig - X_dot_w_reg;
corr_reg_perm = NaN(size(X_dotprod_w_reg_rem, 2), nrep);
for rep_i = 1:nrep
    X_dotprod_w_reg_perm = X_dotprod_w_reg_rem + X_dot_w_reg(perm_idx(:, rep_i), :);
    corr_reg_perm(:, rep_i) = corr(X_dotprod_w_reg_perm, Y_cat);
end
corr_reg_diff = corr_orig - corr_reg_perm;
corr_reg_diff_mean = mean(corr_reg_diff, 2);

%% Results

fimp = struct('nrep', nrep, 'reg', corr_reg_diff_mean);

%% Save feature importance

featimpname = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s.mat', sj_id, parctype));
save(featimpname, 'fimp');

fprintf('Done.\n');

end
