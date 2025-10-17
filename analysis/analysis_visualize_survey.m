function analysis_visualize_survey(varargin)

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

andir = fullfile(bidsdir, 'derivatives', antype);
phenodir = fullfile(bidsdir, 'phenotype');
surdir = fullfile(andir, 'survey');
if ~exist(surdir, 'dir'); mkdir(surdir); end

sur_V1_tab = readtable(fullfile(phenodir, 'survey_V1.tsv'), 'FileType', 'text');

%% Visualize date

hold on;
for sj_num = 1:2
    sj_id = tbl.participant_id(sj_num,:);
    date_mat = table2array(sur_V1_tab(strcmp(sur_V1_tab.participant_id, sj_id), 'timestamp'));
    date_mat = days(date_mat(2:end) - date_mat(2));
    scatter(date_mat, repelem(3-sj_num, 1, numel(date_mat)), 80, 'filled');
end

set(gca, 'FontSize', 14, 'XLim', [-10 255], 'YLim', [0 3], 'XTick', 0:50:250, 'YColor', 'none', ...
    'LineWidth', 2, 'TickLength', [0.01 0.01], 'Tickdir', 'out', 'FontSize', 16);
set(gcf, 'color', 'w');
resize_axes(gca, [850 100]);

figname = fullfile(surdir, 'survey_scandate.pdf');
pagesetup(gcf); saveas(gcf, figname);
close all;

end
