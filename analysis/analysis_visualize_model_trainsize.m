function analysis_visualize_model_trainsize(sj_num, varargin)

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

cols = call_colors;

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat', 'rating_dat_tbin', 'nbin');

trszname = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s.mat', sj_id, parctype));
if ~exist(trszname, 'file'); return; end
load(trszname, 'trsz_nsub', 'trsz_nrep', 'trsz_Yfit', 'trsz_Yfit_tbin');

figbasename = fullfile(preddir, sprintf('%s_predict_rating_trainsize_%s', sj_id, parctype));

%% Time-averaged data

x{1} = reshape(cat(1, rating_dat_tbin{:,:,nbin==10}), [10 size(rating_dat_tbin,1:2)]);
x{2} = reshape(cat(1, rating_dat_tbin{:,:,nbin==5}), [5 size(rating_dat_tbin,1:2)]);
x{3} = cellfun(@mean, rating_dat);
x{4} = mean(x{3});

y{1} = reshape(cat(1, trsz_Yfit_tbin{:,:,:,:,nbin==10}), [10 size(trsz_Yfit_tbin,1:4)]);
y{2} = reshape(cat(1, trsz_Yfit_tbin{:,:,:,:,nbin==5}), [5 size(trsz_Yfit_tbin,1:4)]);
y{3} = squeeze(cellfun(@mean, trsz_Yfit));
y{4} = squeeze(mean(y{3}));

%% Bootstrap samples

rng('default');
nrep = 10000;
boot_idx = cellfun(@(a) NaN(numel(a), nrep), x, 'un', false);

for i = 1:numel(x)
    rep_i = 1;
    while rep_i <= nrep
        boot_idx_each = reshape(1:numel(x{i}), size(x{i}));
        switch i
            case {1 2}
                boot_idx_each = boot_idx_each(:, :, randi(size(boot_idx_each,3), size(boot_idx_each,3), 1));
                for ses_i = 1:size(boot_idx_each,3)
                    boot_idx_each(:,:,ses_i) = boot_idx_each(:, randi(size(boot_idx_each,2), size(boot_idx_each,2), 1), ses_i);
                    for task_i = 1:size(boot_idx_each,2)
                        boot_idx_each(:,task_i,ses_i) = boot_idx_each(randi(size(boot_idx_each,1), size(boot_idx_each,1), 1), task_i, ses_i);
                    end
                end
            case 3
                boot_idx_each = boot_idx_each(:, randi(size(boot_idx_each,2), size(boot_idx_each,2), 1));
                for ses_i = 1:size(boot_idx_each,2)
                    boot_idx_each(:,ses_i) = boot_idx_each(randi(size(boot_idx_each,1), size(boot_idx_each,1), 1), ses_i);
                end
            case 4
                boot_idx_each = boot_idx_each(:, randi(size(boot_idx_each,2), size(boot_idx_each,2), 1));
        end

        boot_idx{i}(:,rep_i) = boot_idx_each(:);
        rep_i = rep_i + 1;
    end
end

%% Correlation all, training size dependency

[rall_r, rall_bootp] = deal(NaN(trsz_nrep*numel(trsz_nsub), numel(x)));
rall_bootci = NaN(trsz_nrep*numel(trsz_nsub), 2, numel(x));

for i = 1:numel(x)
    rall_r(:,i) = corr(x{i}(:), reshape(y{i}, numel(x{i}), []));
    x_boot = x{i}(boot_idx{i});
    y_boot = reshape(y{i}, numel(x{i}), []);
    y_boot = reshape(y_boot(boot_idx{i}, :), [size(x_boot, 1:2) size(y_boot, 2)]);
    rall_boot = squeeze(sum(repmat(zscore(x_boot), [1 1 size(y_boot,3)]) .* zscore(y_boot))).' ./ (size(x_boot,1)-1);
    [rall_bootci(:,:,i), rall_bootp(:,i)] = boot2cip(rall_boot, 'perc', 0);
end

rall_r = reshape(rall_r, [trsz_nrep, numel(trsz_nsub), numel(x)]);
rall_bootci = reshape(rall_bootci, [trsz_nrep, numel(trsz_nsub), 2, numel(x)]);
rall_bootp = reshape(rall_bootp, [trsz_nrep, numel(trsz_nsub), numel(x)]);
rall_bootsig = squeeze(mean(rall_bootp < 0.05 & rall_r > 0));

fprintf('For 4 sessions, mean r = %.2f, %.2f, %.2f, %.2f\n', squeeze(mean(rall_r(:,4,:))));
fprintf('For 4 sessions, sig %% = %d, %d, %d, %d\n', round(rall_bootsig(4,:).*100));
for thr = [50 90 95 97.5 99]
    fprintf('Minimum number of sessions to achieve %.1f%%: %d %d %d %d\n', thr, sum(rall_bootsig < thr/100) + 1);
end

rall_tstats = NaN(3, numel(x));
for i = 1:numel(x)
    [~,~,stats] = glmfit(mean(rall_r(:,:,i)), 1:numel(trsz_nsub));
    rall_tstats(:,i) = [stats.beta(2); stats.t(2); stats.p(2)];
    fprintf('Regression with session, b = %.2f, t = %.2f, P-value %.4f\n', rall_tstats(:,i));
end

hold on;
for i = size(rall_r,3):-1:1
    scatter(trsz_nsub, mean(rall_r(:,:,i)), 50, cols.fourunits(i,:), 'filled');
    errorbar(trsz_nsub, mean(rall_r(:,:,i)), std(rall_r(:,:,i))./trsz_nrep.^0.5.*norminv(0.975), ...
        '-', 'CapSize', 0, 'LineWidth', 2, 'Color', cols.fourunits(i,:));
end
set(gca, 'LineWidth', 2, 'XLim', [trsz_nsub(1)-0.5 trsz_nsub(end)+0.5], 'YLim', extlim(mean(rall_r), 0.1), ...
    'XTick', trsz_nsub, 'TickLength', [0.015 0.015], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);
set(gcf, 'color', 'w');

switch sj_num
    case 1
        resize_axes(gca, [616 350]);
    case 2
        resize_axes(gca, [760 350]);
end

figname = sprintf('%s_corr.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname);
close all;

fprintf('Done.\n');

end
