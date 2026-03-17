function analysis_getsourcedata_model_crosstest(tesj_num, trsj_num, varargin)

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
tesj_id = tbl.participant_id(tesj_num,:);
trsj_id = tbl.participant_id(trsj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
teratdir = fullfile(andir, 'rating', tesj_id);
trpreddir = fullfile(andir, 'predict_rating', trsj_id);

%% Load data

ratname = fullfile(teratdir, sprintf('%s_rating.mat', tesj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat', 'rating_dat_tbin', 'nbin');

crtename = fullfile(trpreddir, sprintf('%s_predict_rating_crosstest_%s_%s.mat', trsj_id, tesj_id, parctype));
if ~exist(crtename, 'file'); return; end
load(crtename, 'crte_Yfit', 'crte_Yfit_tbin');

srcdatbasename = fullfile(trpreddir, sprintf('%s_predict_rating_crosstest_%s_%s_srcdat', trsj_id, tesj_id, parctype));

%% Time-averaged data

x{1} = reshape(cat(1, rating_dat_tbin{:,:,nbin==10}), [10 size(rating_dat_tbin,1:2)]);
x{2} = reshape(cat(1, rating_dat_tbin{:,:,nbin==5}), [5 size(rating_dat_tbin,1:2)]);
x{3} = cellfun(@mean, rating_dat);
x{4} = mean(x{3});

y{1} = reshape(cat(1, crte_Yfit_tbin{:,:,nbin==10}), [10 size(crte_Yfit_tbin,1:2)]);
y{2} = reshape(cat(1, crte_Yfit_tbin{:,:,nbin==5}), [5 size(crte_Yfit_tbin,1:2)]);
y{3} = cellfun(@mean, crte_Yfit);
y{4} = mean(y{3});

xall_hvl = cellfun(@(a) a > quantile(a(:), 0.5), x, 'un', false);

%% Correlation all, scatter

unitnames = ["1min", "2min", "run", "session"];
for i = 1:numel(x)
    T = array2table([x{i}(:) y{i}(:)], 'VariableNames', {'actualpain', 'predictedpain'});
    srcdatname = sprintf('%s_y_yfit_%s.csv', srcdatbasename, unitnames{i});
    writetable(T, srcdatname);
end

%% Classification all, ROC

unitnames = ["1min", "2min", "run", "session"];
for i = 1:numel(x)
    [roc_x, roc_y, ~, ~] = perfcurve(xall_hvl{i}(:), y{i}(:), true);
    T = array2table([roc_x roc_y], 'VariableNames', {'1-specificity', 'sensitivity'});
    srcdatname = sprintf('%s_roc_%s.csv', srcdatbasename, unitnames{i});
    writetable(T, srcdatname);
end

fprintf('Done.\n');

end
