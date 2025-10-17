function analysis_predict_rating_modelingses12(sj_num, varargin)

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
ratdir = fullfile(andir, 'rating', sj_id);
preddir = fullfile(andir, 'predict_rating', sj_id);

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat_rbin', 'nbin');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s.mat', sj_id, parctype));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_rbin');

%% Prepare data

X = func_conn_rbin(:,1:12,nbin==10);
Y = rating_dat_rbin(:,1:12,nbin==10);
idx_train = 1:size(Y,2);
PCs = fliplr(1:numel(cat(1, Y{:,idx_train}))-2);

Yfit_ocv = repmat(Y, 1, 1, numel(PCs));
Yfit_ocv = cellfun(@(a) NaN(size(a)), Yfit_ocv, 'un', false);
Yfit_icv = repmat(Y, 1, 1, numel(PCs), size(Y,2));
Yfit_icv = cellfun(@(a) NaN(size(a)), Yfit_icv, 'un', false);

[beta, beta_edge] = deal([]);
[beta.all, beta_edge.all] = deal(cell(1, numel(PCs)));
[beta.ocv, beta_edge.ocv] = deal(cell(size(Y,2), numel(PCs)));

%% Training: All

wh_tr = idx_train;

fprintf('Train %s ...\n', sprintf('%d,', wh_tr));

X_tr = cat(1, X{:,wh_tr});
Y_tr = cat(1, Y{:,wh_tr});

[func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
X_tr = X_tr*func_PC;
regfit_rocha = lasso_rocha(Y_tr, X_tr(:,1:end-1));
for f = {'beta', 'lambda'}; regfit_rocha.(f{1}) = flipud(regfit_rocha.(f{1})(regfit_rocha.df>0, :)); end
[~, N_rocha] = min(abs(sum(regfit_rocha.beta~=0, 2) - PCs));
L_rocha = regfit_rocha.lambda(N_rocha);
L_matlab = L_rocha .* (numel(Y_tr)/(numel(Y_tr)-1)).^0.5 .* std(Y_tr) + eps;
[regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
assert(isequal(sum(regfit_rocha.beta(N_rocha,:)~=0, 2).', sum(regfit_B~=0)));
regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
beta.all = num2cell(regfit_B, 1);
for p = 1:numel(PCs)
    if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
        beta.all{p} = beta.all{p-1};
    else
        beta.all{p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
    end
end

bcat = cat(2, beta.all{:});
beta_edge.all = num2cell([bcat(1,:); func_PC(:,1:size(bcat,1)-1) * bcat(2:end,:)], 1);

%% Training: OCV

for k = idx_train
    
    wh_tr = setdiff(idx_train, k);
    wh_te = k;
    
    fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));
    
    X_tr = cat(1, X{:,wh_tr});
    X_te = cat(1, X{:,wh_te});
    Y_tr = cat(1, Y{:,wh_tr});

    [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
    X_tr = X_tr*func_PC;
    X_te = X_te*func_PC;
    [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
    regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
    beta.ocv(k,:) = num2cell(regfit_B, 1);
    for p = 1:numel(PCs)
        if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
            beta.ocv{k,p} = beta.ocv{k,p-1};
        else
            beta.ocv{k,p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
        end
    end

    for p = 1:numel(PCs)
        Yfit_ocv(:,wh_te,p) = reshape(...
            mat2cell([ones(size(X_te,1),1), X_te(:,1:numel(beta.ocv{k,p})-1)] * beta.ocv{k,p}, ...
            reshape(cellfun(@numel,Yfit_ocv(:,wh_te,p)),[],1)), ...
            size(Yfit_ocv,1), []);
    end
    bcat = cat(2, beta.ocv{k,:});
    beta_edge.ocv(k,:) = num2cell([bcat(1,:); func_PC(:,1:size(bcat,1)-1) * bcat(2:end,:)], 1);
    
end

%% Training: ICV

for k = idx_train
    
    for kk = setdiff(idx_train, k)
    
        wh_tr = setdiff(idx_train, [k kk]);
        wh_te = kk;
        
        fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));
        
        X_tr = cat(1, X{:,wh_tr});
        X_te = cat(1, X{:,wh_te});
        Y_tr = cat(1, Y{:,wh_tr});
        bb = [];

        [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
        X_tr = X_tr*func_PC;
        X_te = X_te*func_PC;
        [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
        regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
        bb = num2cell(regfit_B, 1);
        for p = 1:numel(PCs)
            if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
                bb{p} = bb{p-1};
            else
                bb{p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
            end
        end
        
        for p = 1:numel(PCs)
            Yfit_icv(:,wh_te,p,k) = reshape(...
                mat2cell([ones(size(X_te,1),1), X_te(:,1:numel(bb{p})-1)] * bb{p}, ...
                reshape(cellfun(@numel,Yfit_icv(:,wh_te,p,k)),[],1)), ...
                size(Yfit_icv,1), []);
        end
        
    end
    
end

%% Calculate results

res = [];

res.L_rocha = L_rocha;
res.L_matlab = L_matlab;

res.ocv_r = NaN(size(Yfit_ocv,3), 1);
for p = 1:numel(PCs)
    res.ocv_r(p) = corr(cat(1, Y{:,:}), cat(1, Yfit_ocv{:,:,p}));
end
[~, res.ocv_bestidx] = max(res.ocv_r);

res.icv_r = NaN(size(Yfit_icv,3:4));
for k = idx_train
    for p = 1:numel(PCs)
        res.icv_r(p,k) = corr(cat(1, Y{:,:}), cat(1, Yfit_icv{:,:,p,k}), 'rows', 'complete');
    end
end
[~, res.icv_bestidx] = max(res.icv_r);

[beta_best, beta_best_edge] = deal([]);
beta_best.all = beta.all{res.ocv_bestidx};
beta_best_edge.all = beta_edge.all{res.ocv_bestidx};
[beta_best.ocv, beta_best_edge.ocv] = deal(cell(size(Y,2), 1));
for k = idx_train
    beta_best.ocv(k,1) = beta.ocv(k,res.icv_bestidx(k));
    beta_best_edge.ocv(k,1) = beta_edge.ocv(k,res.icv_bestidx(k));
end

%% Save model

out = struct('Y', {Y}, 'Yfit_ocv', {Yfit_ocv}, 'Yfit_icv', {Yfit_icv}, ...
    'res', {res}, 'beta', {beta}, 'beta_best', {beta_best}, 'beta_best_edge', {beta_best_edge});
mdlname = fullfile(preddir, sprintf('%s_predict_rating_modelingses12_%s.mat', sj_id, parctype));
save(mdlname, 'out');

fprintf('Done.\n');

end
