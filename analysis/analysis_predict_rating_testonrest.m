function analysis_predict_rating_testonrest(sj_num, varargin)

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

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

tasklabel = 'task-rest';
func = [];
funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabel));
funcmetafile = split(deblank(ls(funcmetafile)));
for ses_i = 1:numel(funcmetafile)
    tempstruct = load(funcmetafile{ses_i}, 'func');
    func{ses_i} = change_bidsdir(tempstruct.func, bidsdir);
end

TR = func{1}.raw.jsondat.RepetitionTime;
ndvol = ceil(30/TR);
yhdel = 13;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

%% Apply model

restYfit = cell(1, numel(func));

for ses_i = 1:numel(func)

    fprintf('Working on ses %d...\n', ses_i);

    wh_keep = logical(importdata(func{ses_i}.den.cen));
    wh_keep = wh_keep(1+ndvol+yhdel:end);

    tempcf = ft_read_cifti_mod(func{ses_i}.cf.den.sbold);
    func_dat = tempcf.data;
    func_dat = func_dat(:, 1+ndvol+yhdel:end).';
    func_dat = func_dat(wh_keep, :);
    func_dat = splitapply(@(x) mean(x,2), func_dat, parc.id.');
    func_dat = func_dat(:, parc.whincl);

    func_conn = zscore(func_dat);
    func_conn = repmat(func_conn, 1, 1, size(func_conn,2));
    func_conn = func_conn .* permute(func_conn, [1 3 2]);
    func_conn = func_conn(:, triu(true(size(func_conn, 2:3)), 1));

    w = out.beta_best_edge.all;
    restYfit{ses_i} = func_conn * w(2:end) + w(1);

end

%% Save results

testonrestname = fullfile(preddir, sprintf('%s_predict_rating_testonrest_%s.mat', sj_id, parctype));
save(testonrestname, 'restYfit');

fprintf('Done.\n');

end