function analysis_parcellation(sj_num, sctxparcname, sctxparclabelname, varargin)

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
            case {'ctxparcname'}
                ctxparcname = varargin{i+1};
            case {'ctxparclabelname'}
                ctxparclabelname = varargin{i+1};
        end
    end
end

addpath(genpath(fullfile(bidsdir, 'code', 'Functions')));

%% Basic setting

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);

andir = fullfile(bidsdir, 'derivatives', antype);
mscdir = fullfile(andir, 'MSC', sj_id);
scandir = fullfile(andir, 'SCAN', sj_id);
parcdir = fullfile(andir, 'parcellation', sj_id);
if ~exist(parcdir, 'dir'); mkdir(parcdir); end

prepdir = fullfile(bidsdir, 'derivatives', preptype);
metadir = fullfile(prepdir, 'metadata');

anatmetafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, 'ses-01', 'anat'));
load(anatmetafile, 'anat');
anat = change_bidsdir(anat, bidsdir);

tasklabel = 'task-rest';
func = [];
funcmetafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, 'ses-*', 'func', tasklabel));
funcmetafile = split(deblank(ls(funcmetafile)));
for ses_i = 1:numel(funcmetafile)
    tempstruct = load(funcmetafile{ses_i}, 'func');
    func{ses_i} = change_bidsdir(tempstruct.func, bidsdir);
end

%% Prepare parcellation

parc.cf = ft_read_cifti_mod(func{1}.cf.den.sbold);
parc.cf.data = [];
parc.cf.dimord = 'scalar_pos';

fsLRdir = fullfile(prepdir, 'ciftify', anat.basename, 'MNINonLinear', 'fsaverage_LR32k');
surf_L = fullfile(fsLRdir, sprintf('%s.L.midthickness.32k_fs_LR.surf.gii', anat.basename));
surf_R = fullfile(fsLRdir, sprintf('%s.R.midthickness.32k_fs_LR.surf.gii', anat.basename));
parc.cf.pos(1:64984,:) = [getfield(gifti(surf_L), 'vertices'); getfield(gifti(surf_R), 'vertices')];

if startsWith(parctype, 'indparc')
    ctxparcname = fullfile(mscdir, sprintf('%s_parcels_edgethresh_0.5.dtseries.nii', sj_id));
    ctxparc_cf = ft_read_cifti_mod(ctxparcname);
    [~, ~, ctxparc_cf.data] = unique(ctxparc_cf.data);
    ctxparc_cf.data = ctxparc_cf.data - 1; % b/c zeros are changed into ones in the unique step above
    ctxparc_cf.data(ctxparc_cf.data == 0) = NaN;
elseif startsWith(parctype, 'schaefer')
    ctxparc_cf = ft_read_cifti_mod(ctxparcname);
    ctxparc_cf.data(ctxparc_cf.data == 0) = NaN;
    ctxparc_cf.data = ctxparc_cf.data(parc.cf.brainstructure(1:numel(ctxparc_cf.data)) ~= -1);
end

tmpdir = tempname(tempdir); mkdir(tmpdir);
system(sprintf('cp %s %s', func{1}.cf.den.sbold, fullfile(tmpdir, 'template.dtseries.nii')));
system(sprintf('wb_command -cifti-create-dense-from-template %s %s -volume-all %s', ...
    fullfile(tmpdir, 'template.dtseries.nii'), fullfile(tmpdir, 'sctxatlas.dscalar.nii'), sctxparcname));
sctxparc_cf = ft_read_cifti_mod(fullfile(tmpdir, 'sctxatlas.dscalar.nii'));
sctxparc_cf.data(sctxparc_cf.data == 0) = NaN;
sctxparc_cf.data = sctxparc_cf.data + max(ctxparc_cf.data);

brainstr_nomw = parc.cf.brainstructure(parc.cf.brainstructure ~= -1);
parc.id = zeros(numel(brainstr_nomw), 1);
parc.id(ismember(brainstr_nomw, 1:2)) = ctxparc_cf.data;
parc.id(~ismember(brainstr_nomw, 1:2)) = sctxparc_cf.data(~ismember(brainstr_nomw, 1:2));

parc.medid = NaN(max(parc.id), 1);
parc.medpos = NaN(max(parc.id), 3);
all_pos = parc.cf.pos(parc.cf.brainstructure ~= -1, :);
for i = 1:max(parc.id)
    roi_idx = find(parc.id == i);
    roi_pos = all_pos(roi_idx, :);
    roi_dist = squeeze(sum((roi_pos.' - permute(roi_pos.', [1 3 2])).^2).^0.5);
    [~, wh_med] = min(sum(roi_dist));
    parc.medidx(i) = roi_idx(wh_med);
    parc.medpos(i,:) = roi_pos(wh_med, :);
end

s = load(sctxparclabelname); f = fieldnames(s);
sctxparclabel = s.(f{1});
sctxparcnet = sctxparclabel.dat(:,2);
sctxparcnetnames = {'THA'; 'HCAMY'; 'BG'; 'CB'; 'BS'};
sctxparcnames = sctxparclabel.names;

if startsWith(parctype, 'indparc')
    ctxparclabelname = fullfile(scandir, sprintf('%s_rawassn_minsize4_regularized_recolored_addSCAN.txt', sj_id));
    ctxparclabel = importdata(ctxparclabelname);
    ctxparcnet = ctxparclabel;
    ctxparcnet(ctxparcnet < 1 | ctxparcnet > 17) = 0;
    ctxparcnetnames = {'NONE'; 'DMN'; 'SCAN'; 'LVIS'; 'FPN'; 'MVIS'; 'DAN'; 'PMOT'; 'LANG'; 'SAL'; 'AMN'; ...
        'HSMN'; 'FSMN'; 'AUD'; 'AMTL'; 'PMTL'; 'PMN'; 'CAN'; 'LSMN'};
    ctxparcnames = strcat({'Ind Parc '}, cellstr(string(1:numel(ctxparcnet)).'));
    parc.netorig = [ctxparcnet; sctxparcnet+17];
    [~, parc.net] = ismember(parc.netorig, [0:1 1.5 2:22]);
elseif startsWith(parctype, 'schaefer')
    ctxparclabel = fileread(ctxparclabelname);
    ctxparclabel = reshape(splitlines(deblank(ctxparclabel)), 2, []).';
    ctxparcnet = regexp(ctxparclabel(:,1), '7Networks_[RL]H_([A-Za-z]+)_', 'tokens');
    [~, ~, ctxparcnet] = unique(cellfun(@(a) a{1}{1}, ctxparcnet, 'un', false), 'stable');
    ctxparcnetnames = {'VN'; 'SMN'; 'DAN'; 'VAN'; 'LN'; 'FPN'; 'DMN';};
    ctxparcnames = ctxparclabel(:,1);
    [parc.netorig, parc.net] = deal([ctxparcnet; sctxparcnet+7]);
end

parc.names = [ctxparcnames; sctxparcnames];
parc.netnames = [ctxparcnetnames; sctxparcnetnames];

if endsWith(parctype, 'wvis')
    parc.whincl = true(size(parc.net));
else
    aparcfile = fullfile(fsLRdir, sprintf('%s.aparc.32k_fs_LR.dlabel.nii', anat.basename));
    aparccf = ft_read_cifti_mod(aparcfile);
    wh_visaparc = ismember((1:numel(parc.net)).', unique(parc.id(find(ismember(aparccf.data, [5 7 11 13 21 40 42 46 48 56])))));
    wh_visnet = ismember(parc.net, find(ismember(parc.netnames, {'LVIS', 'MVIS', 'VN'})));
    if endsWith(parctype, 'onlyvis')
        parc.whincl = wh_visaparc | wh_visnet;
    else
        parc.whincl = ~wh_visaparc & ~wh_visnet;
    end
end

%% Save parcellation

parcname = fullfile(parcdir, sprintf('%s_parcellation_%s_meta.mat', sj_id, parctype));
save(parcname, 'parc');

fprintf('Done.\n');

end
