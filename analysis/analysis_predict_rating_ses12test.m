function analysis_predict_rating_ses12test(sj_num, varargin)

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
ratdir = fullfile(andir, 'rating', sj_id);
preddir = fullfile(andir, 'predict_rating', sj_id);

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
func = func(:,1:12);

TR = func{1,1}.raw.jsondat.RepetitionTime;
ndvol = ceil(30/TR);
yhdel = 13;

%% Load data

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
if ~exist(parcname, 'file'); return; end
load(parcname, 'parc');

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'nbin', 'rbin_idx', 'tbin_idx');

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modelingses12_%s.mat', sj_id, parctype));
if ~exist(mdlname, 'file'); return; end
load(mdlname, 'out');

%% Apply model

Yfit = cell(size(func));
[Yfit_rbin, Yfit_tbin] = deal(cell([size(func) numel(nbin)]));

for task_i = 1:size(func,1)
    
    for ses_i = 1:size(func,2)
        
        fprintf('Working on run %d, ses %d...\n', task_i, ses_i);

        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        wh_keep = wh_keep(1+ndvol+yhdel:end);
        
        tempcf = ft_read_cifti_mod(func{task_i,ses_i}.cf.den.sbold);
        func_dat = tempcf.data;
        func_dat = func_dat(:, 1+ndvol+yhdel:end).';
        func_dat = func_dat(wh_keep, :);
        func_dat = splitapply(@(x) mean(x,2), func_dat, parc.id.');
        func_dat = func_dat(:, parc.whincl);
        
        func_conn = zscore(func_dat);
        func_conn = repmat(func_conn, 1, 1, size(func_conn,2));
        func_conn = func_conn .* permute(func_conn, [1 3 2]);
        func_conn = func_conn(:, triu(true(size(func_conn, 2:3)), 1));

        w = out.beta_best_edge.ocv{ses_i};
        Yfit{task_i, ses_i} = func_conn * w(2:end) + w(1);
        
        for bin_i = 1:numel(nbin)

            Yfit_rbin{task_i,ses_i,bin_i} = splitapply(@(x) mean(x,1), Yfit{task_i, ses_i}, rbin_idx{task_i,ses_i,bin_i});
            Yfit_tbin{task_i,ses_i,bin_i} = splitapply(@(x) mean(x,1), Yfit{task_i, ses_i}, tbin_idx{task_i,ses_i,bin_i});

        end

    end
    
end

%% Save results

testname = fullfile(preddir, sprintf('%s_predict_rating_ses12test_%s.mat', sj_id, parctype));
save(testname, 'Yfit', 'Yfit_rbin', 'Yfit_tbin');

fprintf('Done.\n');

end
