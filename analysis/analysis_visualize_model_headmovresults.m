function analysis_visualize_model_headmovresults(sj_num, opt, varargin)

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
load(obsname, 'nbin');

headmovname = fullfile(preddir, sprintf('%s_predict_rating_headmov.mat', sj_id));
if ~exist(headmovname, 'file'); return; end
load(headmovname, 'fd', 'fd_tbin');

testname = fullfile(preddir, sprintf('%s_predict_rating_test_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
if ~exist(testname, 'file'); return; end
load(testname, 'Yfit', 'Yfit_tbin');

figbasename = strrep(strrep(testname, '_predict_rating_test_', '_predict_rating_figure_headmov_'), '.mat', '');

rcols = call_colors;
ucols = rcols.fourunits;

%% Hold-out

if ~isempty(out.idx_hout)
    Yfit = Yfit(:,out.idx_hout);
    fd = fd(:,out.idx_hout);
    Yfit_tbin = Yfit_tbin(:,out.idx_hout,:);
    fd_tbin = fd_tbin(:,out.idx_hout,:);
end

%% Time-averaged data

x{1} = reshape(cat(1, fd_tbin{:,:,nbin==10}), [10 size(fd_tbin,1:2)]);
x{2} = reshape(cat(1, fd_tbin{:,:,nbin==5}), [5 size(fd_tbin,1:2)]);
x{3} = cellfun(@mean, fd);
x{4} = mean(x{3});

y{1} = reshape(cat(1, Yfit_tbin{:,:,nbin==10}), [10 size(Yfit_tbin,1:2)]);
y{2} = reshape(cat(1, Yfit_tbin{:,:,nbin==5}), [5 size(Yfit_tbin,1:2)]);
y{3} = cellfun(@mean, Yfit);
y{4} = mean(y{3});

xall_hvl = cellfun(@(a) a > quantile(a(:), 0.5), x, 'un', false);

s = calc_res(squeeze(num2cell(x{1},1)), squeeze(num2cell(y{1},1)));

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

        if numel(unique(xall_hvl{i}(boot_idx_each))) == 1; disp('Skip boot iter'); continue; end
        boot_idx{i}(:,rep_i) = boot_idx_each(:);
        rep_i = rep_i + 1;
    end
end

x_boot = cellfun(@(a,b) squeeze(reshape(a(b), [size(a) size(b,2)])), x, boot_idx, 'un', false);
y_boot = cellfun(@(a,b) squeeze(reshape(a(b), [size(a) size(b,2)])), y, boot_idx, 'un', false);

%% Correlation all, scatter

[rall_r, rall_bootp] = deal(NaN(1, numel(x)));
rall_bootci = NaN(2, numel(x));

unitnames = {'1min', '2min', 'run', 'ses'};
sc_sz = [40 20 20 20];
for i = 1:numel(x)
    rall_r(i) = corr(x{i}(:), y{i}(:));
    rall_boot = cellfun(@corr, ...
        num2cell(reshape(x_boot{i},[],nrep),1), num2cell(reshape(y_boot{i},[],nrep),1));
    [rall_bootci(:,i), rall_bootp(i)] = boot2cip(rall_boot, 'perc', 0);
    fprintf('r = %.2f (95%% CI %.2f-%.2f), P-value %.4f\n', rall_r(i), rall_bootci(:,i), rall_bootp(i));
    
    hold on;
    sc_a = scatter(x{i}(:), y{i}(:), sc_sz(i), ucols(i,:), 'filled', 'MarkerFaceAlpha', 0.8);
    set(gca, 'LineWidth', 2, 'TickLength', [0.03 0.03], 'Tickdir', 'out', 'FontSize', 16);
    switch sj_num
        case 1
            switch i
                case 1
                    set(gca, 'XLim', [0 0.26], 'YLim', [0.4 0.8], 'XTick', 0:0.05:0.25, 'YTick', 0.4:0.1:0.8);
                case 2
                    set(gca, 'XLim', [0.03 0.21], 'YLim', [0.5 0.75], 'XTick', 0.06:0.06:0.18, 'YTick', 0.5:0.1:0.7);
                case 3
                    set(gca, 'XLim', [0.05 0.19], 'YLim', [0.5 0.75], 'XTick', 0.06:0.06:0.18, 'YTick', 0.5:0.1:0.7);
                case 4
                    set(gca, 'XLim', [0.05 0.19], 'YLim', [0.5 0.75], 'XTick', 0.06:0.06:0.18, 'YTick', 0.5:0.1:0.7);
            end
        case 2
            switch i
                case 1
                    set(gca, 'XLim', [0 0.3], 'YLim', [0 0.8], 'XTick', 0:0.05:0.3, 'YTick', 0:0.2:0.8);
                case 2
                    set(gca, 'XLim', [0 0.3], 'YLim', [0 0.8], 'XTick', 0:0.1:0.3, 'YTick', 0:0.4:0.8);
                case {3 4}
                    set(gca, 'XLim', [0.07 0.21], 'YLim', [0 0.8], 'XTick', 0.1:0.05:0.2, 'YTick', 0:0.4:0.8);
            end
    end
    scref_a = line(xlim, glmfit(x{i}(:), y{i}(:)).'*[ones(1,2); xlim], 'LineWidth', 4, 'Color', [0.5 0.5 0.5 1]);
    set(gca, 'Children', [sc_a scref_a]);

    switch i
        case 1
            resize_axes(gca, [300 300]);
        case {2 3 4}
            resize_axes(gca, [120 120]);
    end
    set(gcf, 'color', 'w');

    figname = sprintf('%s_y_yfit_%s_scatter.pdf', figbasename, unitnames{i});
    pagesetup(gcf); saveas(gcf, figname); saveas(gcf, strrep(figname, '.pdf', '.png'));
    close all;
end

%% Classification all, ROC

[rocall_auc, rocall_bootci, rocall_bootp] = deal(NaN(1,4), NaN(2,4), NaN(1,4));

roc_a = gobjects(1,numel(x));
hold on;
plot([0 1], [0 1], 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
for i = 1:numel(x)
    [roc_x, roc_y, ~, rocall_auc(i)] = perfcurve(xall_hvl{i}(:), y{i}(:), true);
    roc_a(i) = plot(roc_x, roc_y, 'LineWidth', 2, 'Color', ucols(i,:));
    [~, ~, ~, rocall_boot] = cellfun(@(a,b) perfcurve(a, b, true), ...
        num2cell(xall_hvl{i}(boot_idx{i}),1), num2cell(reshape(y_boot{i},[],nrep),1), 'un', false);
    [rocall_bootci(:,i), rocall_bootp(i)] = boot2cip(cell2mat(rocall_boot), 'perc', 0.5);
    fprintf('AUC = %.2f (95%% CI %.2f-%.2f), P-value %.4f\n', rocall_auc(i), rocall_bootci(:,i), rocall_bootp(i));
end

set(gca, 'LineWidth', 2, 'XLim', [-0.005 1], 'YLim', [-0.005 1], ...
    'TickLength', [0.03 0.03], 'Tickdir', 'out', 'Box', 'off', 'FontSize', 16);

resize_axes(gca, [300 300]);
set(gcf, 'color', 'w');

figname = sprintf('%s_roc_all_plot.pdf', figbasename);
pagesetup(gcf); saveas(gcf, figname); saveas(gcf, strrep(figname, '.pdf', '.png'));
close all;

end