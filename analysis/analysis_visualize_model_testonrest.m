function analysis_visualize_model_testonrest(sj_num, varargin)

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
load(ratname, 'rating_dat');

testonrestname = fullfile(preddir, sprintf('%s_predict_rating_testonrest_%s.mat', sj_id, parctype));
if ~exist(testonrestname, 'file'); return; end
load(testonrestname, 'restYfit');

figbasename = fullfile(preddir, sprintf('%s_predict_rating_testonrest_%s', sj_id, parctype));

%% Time-averaged data

x = cellfun(@mean, rating_dat(1,:));
y = cellfun(@mean, restYfit);

xall_hvl = x > quantile(x(:), 0.5);

%% Bootstrap samples

rng('default');
nrep = 10000;
boot_idx = randi(numel(x), numel(x), nrep);
x_boot = x(boot_idx);
y_boot = y(boot_idx);

%% Correlation all, scatter

rall_r = corr(x(:), y(:));
rall_boot = cellfun(@corr, num2cell(x_boot,1), num2cell(y_boot,1));
[rall_bootci, rall_bootp] = boot2cip(rall_boot, 'perc', 0);
fprintf('r = %.2f (95%% CI %.2f-%.2f), P-value %.4f\n', rall_r, rall_bootci, rall_bootp);

hold on;
sc_a = scatter(x(:), y(:), 40, cols.fourunits(4,:), 'filled', 'MarkerFaceAlpha', 0.8);
set(gca, 'LineWidth', 2, 'TickLength', [0.03 0.03], 'Tickdir', 'out', 'FontSize', 16);
switch sj_num
    case 1
        set(gca, 'XLim', [0.3 0.75], 'YLim', [0.35 0.6], 'XTick', 0.3:0.1:0.7, 'YTick', 0.4:0.1:0.6);
    case 2
        set(gca, 'XLim', [0 0.67], 'YLim', [0 0.67], 'XTick', 0:0.2:0.6, 'YTick', 0:0.2:0.6);
end
scref_a = line(xlim, glmfit(x(:), y(:)).'*[ones(1,2); xlim], 'LineWidth', 4, 'Color', [0.5 0.5 0.5 1]);
set(gca, 'Children', [sc_a scref_a]);

resize_axes(gca, [300 300]);
set(gcf, 'color', 'w');

figname = sprintf('%s_y_yfit_scatter.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname);
close all;

%% Classification all, ROC

hold on;
plot([0 1], [0 1], 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);

[roc_x, roc_y, ~, rocall_auc] = perfcurve(xall_hvl(:), y(:), true);
plot(roc_x, roc_y, 'LineWidth', 2, 'Color', cols.fourunits(4,:));
[~, ~, ~, rocall_boot] = cellfun(@(a,b) perfcurve(a, b, true), ...
    num2cell(xall_hvl(boot_idx),1), num2cell(reshape(y_boot,[],nrep),1), 'un', false);
[rocall_bootci, rocall_bootp] = boot2cip(cell2mat(rocall_boot), 'perc', 0.5);
fprintf('AUC = %.2f (95%% CI %.2f-%.2f), P-value %.4f\n', rocall_auc, rocall_bootci, rocall_bootp);

set(gca, 'LineWidth', 2, 'XLim', [-0.005 1], 'YLim', [-0.005 1], ...
    'TickLength', [0.03 0.03], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);

resize_axes(gca, [300 300]);
set(gcf, 'color', 'w');

figname = sprintf('%s_roc_all_plot.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname);
close all;

fprintf('Done.\n');

end
