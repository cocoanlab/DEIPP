function analysis_predict_rating_func_featimp(sj_num, opt, varargin)

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

%% Load data

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

obsname = fullfile(preddir, sprintf('%s_predict_rating_observation.mat', sj_id));
if ~exist(obsname, 'file'); return; end
load(obsname, 'rating_dat_tbin', 'nbin');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s_%s.mat', ...
    sj_id, opt.dfc, opt.parc));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_tbin');

%%

X = func_conn_tbin(:,:,nbin==opt.nbin);
Y = rating_dat_tbin(:,:,nbin==opt.nbin);

w = repmat({out.beta_best_edge.all}, 1, size(X,2));
w(out.idx_train) = out.beta_best_edge.ocv;

if ~isempty(out.idx_hout)
    X = X(:,out.idx_hout);
    Y = Y(:,out.idx_hout);
    w = w(out.idx_hout);
end

Y_cat = cat(1, Y{:});

w_int = cellfun(@(a) a(1, :), w);
w_beta = cellfun(@(a) a(2:end, :), w, 'un', false);

X_dot_w = cellfun(@(a,b) a .* b.', X, repmat(w_beta, size(X,1), 1), 'un', false);
X_dot_w = cat(1, X_dot_w{:});
X_prod_w_orig = sum(X_dot_w, 2) + reshape(repmat(w_int, opt.nbin*size(X,1), 1), [], 1);
corr_orig = corr(X_prod_w_orig, Y_cat);
mse_orig = mean((X_prod_w_orig - Y_cat).^2);

rng('default');
nrep = 10000;
perm_idx = NaN(numel(Y_cat), nrep);
% for rep_i = 1:nrep
%     perm_idx(:,rep_i) = randperm(size(perm_idx,1));
% end
for rep_i = 1:nrep
    perm_idx_each = reshape(1:size(perm_idx,1), opt.nbin, size(Y,1), size(Y,2));
    perm_idx_each = perm_idx_each(:, :, randperm(size(perm_idx_each,3)));
    for ses_i = 1:size(perm_idx_each,3)
        perm_idx_each(:,:,ses_i) = perm_idx_each(:, randperm(size(perm_idx_each,2)), ses_i);
        for task_i = 1:size(perm_idx_each,2)
            perm_idx_each(:,task_i,ses_i) = perm_idx_each(randperm(size(perm_idx_each,1)), task_i, ses_i);
        end
    end
    perm_idx(:,rep_i) = perm_idx_each(:);
end

% Region
X_dot_w_reg = cellfun(@refmt_r, num2cell(X_dot_w, 2), 'un', false);
X_dot_w_reg = cellfun(@sum, X_dot_w_reg, 'un', false);
X_dot_w_reg = cat(1, X_dot_w_reg{:});
X_prod_w_reg_rem = X_prod_w_orig - X_dot_w_reg;
[corr_reg_perm, mse_reg_perm] = deal(NaN(size(X_prod_w_reg_rem, 2), nrep));
fprintf('\nWorking on feature importance:\n');
for rep_i = 1:nrep
    if mod(rep_i,1000)==0; fprintf('Working on feature importance iter %.6d / %.6d.\n', rep_i, nrep); end
    X_prod_w_reg_perm = X_prod_w_reg_rem + X_dot_w_reg(perm_idx(:,rep_i), :);
    corr_reg_perm(:,rep_i) = corr(X_prod_w_reg_perm, Y_cat);
    mse_reg_perm(:,rep_i) = mean((X_prod_w_reg_perm - Y_cat).^2);
end
corr_reg_diff_all = corr_orig - corr_reg_perm;
mse_reg_diff_all = -(mse_orig - mse_reg_perm);

% Results
fimp.nrep = nrep;
fimp.orig.corr = corr_orig;
fimp.orig.mse = mse_orig;
fimp.reg.corr.all = corr_reg_diff_all;
fimp.reg.corr.mean = mean(corr_reg_diff_all, 2);
fimp.reg.corr.std = std(corr_reg_diff_all, [], 2);
fimp.reg.corr.p_ot = (sum(corr_reg_diff_all <= 0, 2) + 1) ./ (fimp.nrep+1);
fimp.reg.corr.p_tt = min(fimp.reg.corr.p_ot, 1 - fimp.reg.corr.p_ot) .* 2;
fimp.reg.mse.all = mse_reg_diff_all;
fimp.reg.mse.mean = mean(mse_reg_diff_all, 2);
fimp.reg.mse.std = std(mse_reg_diff_all, [], 2);
fimp.reg.mse.p_ot = (sum(mse_reg_diff_all <= 0, 2) + 1) ./ (fimp.nrep+1);
fimp.reg.mse.p_tt = min(fimp.reg.mse.p_ot, 1 - fimp.reg.mse.p_ot) .* 2;

%%

featimpname = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
save(featimpname, 'fimp');

fprintf('Done.\n');

end
