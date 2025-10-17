function analysis_univariate_rating_firstlevel(sj_num, spmdir, varargin)

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

addpath(genpath(spmdir));
addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
ratdir = fullfile(andir, 'rating', sj_id);
univdir = fullfile(andir, 'univariate_rating', sj_id);
if ~exist(univdir, 'dir'); mkdir(univdir); end

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

funcrefcf = ft_read_cifti_mod(func{1,1}.cf.den.sbold);
funcrefcf.data = [];
funcrefcf.dimord = 'scalar_pos';

TR = func{1,1}.raw.jsondat.RepetitionTime;
ndvol = ceil(30/TR);

%% Load data

ratname = fullfile(ratdir, sprintf('%s_rating.mat', sj_id));
if ~exist(ratname, 'file'); return; end
load(ratname, 'rating_dat_orig');

%% First-level

b_all = cell(1, size(func,2));

for ses_i = 1:size(func,2)

    fprintf('Working on ses %d...\n', ses_i);

    [X_e, X_n] = deal(cell(1, size(func,1)));
    Y = [];
    
    for task_i = 1:size(func,1)
        
        rating_dat_hrfconv = conv(rating_dat_orig{task_i,ses_i}, spm_hrf(TR));
        X_e{task_i} = rating_dat_hrfconv(1+ndvol:numel(rating_dat_orig{task_i,ses_i}));
        
        wh_keep = logical(importdata(func{task_i,ses_i}.den.cen));
        ortmat = importdata(func{task_i,ses_i}.den.ort); % intercept, linear trend, HPF 0.005Hz, 6 RP, 10 aCompCor
        X_n{task_i} = [ortmat (find(~wh_keep)==1:size(ortmat,1)).'];
        X_n{task_i} = X_n{task_i}(1+ndvol:end, :);
        X_n{task_i} = X_n{task_i}(:, any(logical(X_n{task_i})));

        tempcf = ft_read_cifti_mod(func{task_i,ses_i}.cf.mni.sbold);
        Y = [Y; tempcf.data(:, 1+ndvol:end).'];
        
    end

    X = [cat(1, X_e{:}) blkdiag(X_n{:})];
    
    b = X \ Y;
    b_all{ses_i} = b(1,:).';

end

b_all = cat(2, b_all{:});

%% Save results

tempcf = funcrefcf;
tempcf.data = b_all;
tempcf.mapname = strcat({'beta, session '}, cellstr(string(1:size(b_all,2)).'));
imgname = fullfile(univdir, sprintf('%s_univariate_rating_firstlevel_b', sj_id));
ft_write_cifti_mod(imgname, tempcf);

fprintf('Done.\n');

end
