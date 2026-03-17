function analysis_getsourcedata_model_featimp(sj_num, varargin)

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
parcdir = fullfile(andir, 'parcellation', sj_id);
preddir = fullfile(andir, 'predict_rating', sj_id);

cols = call_colors;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

featimpname = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s.mat', sj_id, parctype));
if ~exist(featimpname, 'file'); return; end
load(featimpname, 'fimp');

srcdatbasename = fullfile(preddir, sprintf('%s_predict_rating_featimp_%s_srcdat', sj_id, parctype));

%% Prepare feature importance table

fimp_tbl = table((1:numel(parc.net)).', string(parc.names), round(parc.medpos), string(parc.netnames(parc.net)), ...
    'VariableNames', {'ROI', 'name', 'MNI', 'netnames'});
fimp_tbl = fimp_tbl(parc.whincl, :);
fimp_tbl = addvars(fimp_tbl, fimp.reg, 'After', 1, 'NewVariableNames', 'imp');

srcdatname = sprintf('%s.csv', srcdatbasename);
writetable(fimp_tbl, srcdatname);

fprintf('Done.\n');

end
