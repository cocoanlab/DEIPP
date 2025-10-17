function analysis_visualize_rating(sj_num, varargin)

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
ratdir = fullfile(andir, 'rating', sj_id);

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
nvol = func{1,1}.rmd.nvol;
ndvol = ceil(30/TR);

cols = call_colors;

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat_orig', 'rating_dat', 'rating_dat_tbin', 'nbin');

figbasename = fullfile(ratdir, sprintf('%s_rating', sj_id));

%% TR-level rating

for ses_i = 1:size(func,2)
    for task_i = 1:size(func,1)
        plot((0:(nvol-ndvol-1))*TR/60, rating_dat_orig{task_i,ses_i}(1+ndvol:end), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        set(gca, 'XLim', [-1 11], 'YLim', [-0.075 1.075], 'XTick', 0:2:10, 'YTick', 0:0.2:1, ...
            'TickDir', 'out', 'TickLength', [0.015 0.015], 'XColor', 'none', 'YColor', 'none', ...
            'Box', 'off', 'XGrid', 'on', 'YGrid', 'on', 'GridLineWidth', 1.5);

        resize_axes(gca, [50 75]);
        set(gcf, 'color', 'w');

        figname = sprintf('%s_TR_ses-%.2d_run-%d.pdf', figbasename, ses_i, task_i);
        exportgraphics(gcf, figname, 'BackgroundColor', 'none');
        close all;
    end
end

%% Bin-level rating

x{1} = reshape(cat(1, rating_dat_tbin{:,:,nbin==10}), [10 size(rating_dat_tbin,1:2)]);
x{2} = reshape(cat(1, rating_dat_tbin{:,:,nbin==5}), [5 size(rating_dat_tbin,1:2)]);
x{3} = cellfun(@mean, rating_dat);
x{4} = mean(x{3});

tbin_fac = numel(x{1}) ./ cellfun(@numel, x);
unitnames = {'1min', '2min', 'run', 'ses'};
sc_sz = [10 20 40 80];
for i = 1:numel(x)
    fprintf('Mean: %.2f, SD: %.2f\n', mean(x{i}(:)), std(x{i}(:)));
    hold on;
    scatter((1+tbin_fac(i))/2 : tbin_fac(i) : numel(x{1}), x{i}(:), sc_sz(i), cols.fourunits(i,:), 'filled');
    set(gca, 'XLim', [-2 numel(x{1})+3], 'YLim', [0 1], 'XTick', 0.5 : tbin_fac(4) : 0.5+numel(x{1}), 'XTickLabel', [], 'YTick', 0:0.2:1, ...
        'TickDir', 'out', 'TickLength', [0.21 0.21] ./ numel(x{4}), 'Box', 'off', 'FontSize', 16, ...
        'YGrid', 'on', 'GridLineStyle', '--', 'XMinorTick', 'on');
    set(get(gca, 'XAxis'), 'MinorTickValues', 0.5 : tbin_fac(3) : 0.5+numel(x{1}));
    text((1+tbin_fac(4))/2 : tbin_fac(4) : numel(x{1}), repmat(-0.15, 1, numel(x{4})), ...
         split(deblank(sprintf('ses-%.2d\n', 1:numel(x{4})))), 'FontSize', 16, 'HorizontalAlignment', 'center');

    switch sj_num
        case 1
            set(gca, 'Position', [0.03 0.2 0.965 0.75]);
            resize_axes(gca, [1498 100]);
        case 2
            set(gca, 'Position', [0.025 0.2 0.97 0.75]);
            resize_axes(gca, [1820 100]);
        case 3
            set(gca, 'Position', [0.06 0.2 0.935 0.75]);
            resize_axes(gca, [785 100]);
    end
    set(gcf, 'color', 'w');

    figname = sprintf('%s_tbin_%s.pdf', figbasename, unitnames{i});
    pagesetup(gcf); saveas(gcf, figname);
    close all;

    histogram(x{i}(:), 'BinWidth', 0.05, 'FaceColor', cols.fourunits(i,:));
    set(gca, 'XLim', [0 1], 'XTick', 0:0.1:1, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off', 'FontSize', 16);

    resize_axes(gca, [360 90]);
    set(gcf, 'color', 'w');

    figname = sprintf('%s_tbinhist_%s.pdf', figbasename, unitnames{i});
    pagesetup(gcf); saveas(gcf, figname);
    close all;
end

fprintf('Done.\n');

end
