function analysis_getsourcedata_model_trainsize(sj_num, varargin)

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
load(ratname, 'rating_dat', 'rating_dat_tbin', 'nbin');

trszname = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s.mat', sj_id, parctype));
if ~exist(trszname, 'file'); return; end
load(trszname, 'trsz_nsub', 'trsz_nrep', 'trsz_Yfit', 'trsz_Yfit_tbin');

srcdatbasename = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s_srcdat', sj_id, parctype));

%% Time-averaged data

x{1} = reshape(cat(1, rating_dat_tbin{:,:,nbin==10}), [10 size(rating_dat_tbin,1:2)]);
x{2} = reshape(cat(1, rating_dat_tbin{:,:,nbin==5}), [5 size(rating_dat_tbin,1:2)]);
x{3} = cellfun(@mean, rating_dat);
x{4} = mean(x{3});

y{1} = reshape(cat(1, trsz_Yfit_tbin{:,:,:,:,nbin==10}), [10 size(trsz_Yfit_tbin,1:4)]);
y{2} = reshape(cat(1, trsz_Yfit_tbin{:,:,:,:,nbin==5}), [5 size(trsz_Yfit_tbin,1:4)]);
y{3} = squeeze(cellfun(@mean, trsz_Yfit));
y{4} = squeeze(mean(y{3}));

%% Correlation all, training size dependency

rall_r = NaN(trsz_nrep*numel(trsz_nsub), numel(x));

for i = 1:numel(x)
    rall_r(:,i) = corr(x{i}(:), reshape(y{i}, numel(x{i}), []));
end

rall_r = reshape(rall_r, [trsz_nrep, numel(trsz_nsub), numel(x)]);

T = table((1:numel(trsz_nsub)).', 'VariableNames', {'n_session'});
unitnames = ["1min", "2min", "run", "session"];
for i = 1:numel(x)
    T = [T array2table(mean(rall_r(:,:,i)).' + [0 -1 1] .* std(rall_r(:,:,i)).'./trsz_nrep.^0.5.*norminv(0.975), ...
        'VariableNames', unitnames{i} + "_" + ["mean", "95ci_lower", "95ci_upper"])];
end
srcdatname = sprintf('%s_corr.csv', srcdatbasename);
writetable(T, srcdatname);

fprintf('Done.\n');

end
