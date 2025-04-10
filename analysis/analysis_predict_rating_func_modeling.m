function analysis_predict_rating_func_modeling(sj_num, opt, varargin)

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

obsname = fullfile(preddir, sprintf('%s_predict_rating_observation.mat', sj_id));
if ~exist(obsname, 'file'); return; end
load(obsname, 'rating_dat_rbin', 'rating_dat_tbin', 'nbin');

featname = fullfile(preddir, sprintf('%s_predict_rating_feature_%s_%s.mat', ...
    sj_id, opt.dfc, opt.parc));
if ~exist(featname, 'file'); return; end
load(featname, 'func_conn_rbin', 'func_conn_tbin');

%% Prepare data

switch opt.bintype
    case 'rbin'
        X = func_conn_rbin(:,:,nbin==opt.nbin);
        Y = rating_dat_rbin(:,:,nbin==opt.nbin);
    case 'tbin'
        X = func_conn_tbin(:,:,nbin==opt.nbin);
        Y = rating_dat_tbin(:,:,nbin==opt.nbin);
end

switch opt.hout
    case 'none'
        idx_hout = [];
    case 'rand'
        rng('default');
        while true
            idx_hout = unique(randperm(size(Y,2), size(Y,2)-opt.ntrain));
            h = kstest2(cat(1, Y{:,setdiff(1:size(Y,2),idx_hout)}), cat(1, Y{:,idx_hout}));
            if h == 0; break; end
        end
    case 'pro'
        idx_hout = opt.ntrain+1:size(Y,2);
    otherwise
        error('Wrong hout!');
end
idx_train = setdiff(1:size(Y,2), idx_hout);

Yfit_ocv = repmat(Y, 1, 1, numel(opt.param));
Yfit_ocv = cellfun(@(a) NaN(size(a)), Yfit_ocv, 'un', false);
Yfit_icv = repmat(Y, 1, 1, numel(opt.param), size(Y,2));
Yfit_icv = cellfun(@(a) NaN(size(a)), Yfit_icv, 'un', false);

[beta, beta_edge] = deal([]);
[beta.all, beta_edge.all] = deal(cell(1, numel(opt.param)));
[beta.ocv, beta_edge.ocv] = deal(cell(size(Y,2), numel(opt.param)));

%% Training: All

wh_tr = idx_train;
wh_te = idx_hout;

fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));

X_tr = cat(1, X{:,wh_tr});
X_te = cat(1, X{:,wh_te});
Y_tr = cat(1, Y{:,wh_tr});

switch opt.algorithm
    case {'LASSOPCR', 'LASSOOLSPCR'}
        assert(isequal(sort(opt.param, 'descend'), opt.param));
        [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
        X_tr = X_tr*func_PC;
        if ~isempty(wh_te); X_te = X_te*func_PC; end
        regfit_rocha = lasso_rocha(Y_tr, X_tr(:,1:end-1));
        for f = {'beta', 'lambda'}; regfit_rocha.(f{1}) = flipud(regfit_rocha.(f{1})(regfit_rocha.df>0, :)); end
        [~, N_rocha] = min(abs(sum(regfit_rocha.beta~=0, 2) - opt.param));
        L_rocha = regfit_rocha.lambda(N_rocha);
        L_matlab = L_rocha .* (numel(Y_tr)/(numel(Y_tr)-1)).^0.5 .* std(Y_tr) + eps;
        [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
        assert(isequal(sum(regfit_rocha.beta(N_rocha,:)~=0, 2).', sum(regfit_B~=0)));
        regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
        beta.all = num2cell(regfit_B, 1);
        if strcmp(opt.algorithm, 'LASSOOLSPCR')
            for p = 1:numel(opt.param)
                if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
                    beta.all{p} = beta.all{p-1};
                else
                    beta.all{p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
                end
            end
        end
end

for p = 1:numel(opt.param)
    if ~isempty(wh_te)
        Yfit_ocv(:,wh_te,p) = reshape(...
            mat2cell([ones(size(X_te,1),1), X_te(:,1:numel(beta.all{p})-1)] * beta.all{p}, ...
            reshape(cellfun(@numel,Yfit_ocv(:,wh_te,p)),[],1)), ...
            size(Yfit_ocv,1), []);
    end
end
bcat = cat(2, beta.all{:});
beta_edge.all = num2cell([bcat(1,:); func_PC(:,1:size(bcat,1)-1) * bcat(2:end,:)], 1);

%% Training: CV

for k = idx_train
    
    wh_tr = setdiff(idx_train, k);
    wh_te = k;
    
    fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));
    
    X_tr = cat(1, X{:,wh_tr});
    X_te = cat(1, X{:,wh_te});
    Y_tr = cat(1, Y{:,wh_tr});
    
    switch opt.algorithm
        case {'LASSOPCR', 'LASSOOLSPCR'}
            [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
            X_tr = X_tr*func_PC;
            X_te = X_te*func_PC;
            [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
            regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
            beta.ocv(k,:) = num2cell(regfit_B, 1);
            if strcmp(opt.algorithm, 'LASSOOLSPCR')
                for p = 1:numel(opt.param)
                    if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
                        beta.ocv{k,p} = beta.ocv{k,p-1};
                    else
                        beta.ocv{k,p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
                    end
                end
            end
    end
    
    for p = 1:numel(opt.param)
        Yfit_ocv(:,wh_te,p) = reshape(...
            mat2cell([ones(size(X_te,1),1), X_te(:,1:numel(beta.ocv{k,p})-1)] * beta.ocv{k,p}, ...
            reshape(cellfun(@numel,Yfit_ocv(:,wh_te,p)),[],1)), ...
            size(Yfit_ocv,1), []);
    end
    bcat = cat(2, beta.ocv{k,:});
    beta_edge.ocv(k,:) = num2cell([bcat(1,:); func_PC(:,1:size(bcat,1)-1) * bcat(2:end,:)], 1);
    
end

%% Training: NCV

for k = idx_train
    
    for kk = setdiff(idx_train, k)
    
        wh_tr = setdiff(idx_train, [k kk]);
        wh_te = kk;
        
        fprintf('Train %s Test %s ...\n', sprintf('%d,', wh_tr), sprintf('%d,', wh_te));
        
        X_tr = cat(1, X{:,wh_tr});
        X_te = cat(1, X{:,wh_te});
        Y_tr = cat(1, Y{:,wh_tr});
        bb = [];
        
        switch opt.algorithm
            case {'LASSOPCR', 'LASSOOLSPCR'}
                [func_PC, ~, ~] = svd((X_tr - mean(X_tr)).', 'econ');
                X_tr = X_tr*func_PC;
                X_te = X_te*func_PC;
                [regfit_B, regfit_stats] = lasso(X_tr(:,1:end-1), Y_tr, 'Alpha', 1, 'Lambda', L_matlab);
                regfit_B = [regfit_stats.Intercept; regfit_B; zeros(1,size(regfit_B,2))];
                bb = num2cell(regfit_B, 1);
                if strcmp(opt.algorithm, 'LASSOOLSPCR')
                    for p = 1:numel(opt.param)
                        if p > 1 && isequal(regfit_B(:,p)~=0, regfit_B(:,p-1)~=0)
                            bb{p} = bb{p-1};
                        else
                            bb{p}(regfit_B(:,p)~=0) = [ones(size(X_tr,1),1), X_tr(:,regfit_B(2:end,p)~=0)] \ Y_tr;
                        end
                    end
                end
        end
        
        for p = 1:numel(opt.param)
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

res.icv = repmat(calc_res(Y, cellfun(@(a) NaN(size(a)), Y, 'un', false), idx_hout), numel(opt.param), size(Yfit_icv, 4));
for k = 1:size(Yfit_icv, 4)
    for p = 1:numel(opt.param)
        res.icv(p,k) = calc_res(Y, Yfit_icv(:,:,p,k), idx_hout);
    end
end
[~, res.icv_bestidx] = max(reshape(cellfun(@(a) a.('r_r'), {res.icv.train}), size(res.icv)), [], 1);
% [~, res.icv_bestidx] = min(reshape(cellfun(@(a) a.('mse'), {res.icv.train}), size(res.icv)), [], 1);
res.icv_bestidx(idx_hout) = NaN;

res.ocv = repmat(calc_res(Y, cellfun(@(a) NaN(size(a)), Y, 'un', false), idx_hout), numel(opt.param), 1);
for p = 1:numel(opt.param)
    res.ocv(p) = calc_res(Y, Yfit_ocv(:,:,p), idx_hout);
end
[~, res.ocv_bestidx] = max(reshape(cellfun(@(a) a.('r_r'), {res.ocv.train}), size(res.ocv)), [], 1);
% [~, res.ocv_bestidx] = min(reshape(cellfun(@(a) a.('mse'), {res.ocv.train}), size(res.ocv)), [], 1);

Yfit_ncv = cell(size(Y));
for j = 1:size(Y, 2)
    if ismember(j, idx_train)
        Yfit_ncv(:,j) = Yfit_ocv(:,j,res.icv_bestidx(j));
    else
        Yfit_ncv(:,j) = Yfit_ocv(:,j,res.ocv_bestidx);
    end
end

res.ncv = calc_res(Y, Yfit_ncv, idx_hout);

[beta_best, beta_best_edge] = deal([]);
beta_best.all = beta.all{res.ocv_bestidx};
beta_best_edge.all = beta_edge.all{res.ocv_bestidx};
[beta_best.ocv, beta_best_edge.ocv] = deal(cell(size(Y,2), 1));
for j = idx_train
    beta_best.ocv(j,1) = beta.ocv(j,res.icv_bestidx(j));
    beta_best_edge.ocv(j,1) = beta_edge.ocv(j,res.icv_bestidx(j));
end

%% Save model

out = struct('opt', {opt}, 'Y', {Y}, 'Yfit_ncv', {Yfit_ncv}, 'Yfit_ocv', {Yfit_ocv}, 'Yfit_icv', {Yfit_icv}, ...
    'res', {res}, 'beta', {beta}, 'beta_best', {beta_best}, 'beta_best_edge', {beta_best_edge}, 'idx_train', idx_train, 'idx_hout', idx_hout);
mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
save(mdlname, 'out');

fprintf('Done.\n');

end