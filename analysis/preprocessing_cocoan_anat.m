function preprocessing_cocoan_anat(sj_num, ses_num, cfgitdir, oasisdir, varargin)

do_override = false;
preptype = 'cocoan-preproc';
bidsdir = fileparts(fileparts(mfilename('fullpath'))); % mfilename: bidsdir/code/~.m

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'do_override'}
                do_override = true;
            case {'preptype'}
                preptype = varargin{i+1};
            case {'bidsdir'}
                bidsdir = varargin{i+1};
        end
    end
end

%% Basic setting

prepdir = fullfile(bidsdir, 'derivatives', preptype);

tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));

sj_id = tbl.participant_id(sj_num,:);
ses_id = sprintf('ses-%.2d', ses_num);
sesnext_id = sprintf('ses-%.2d', ses_num+1); % for T2

logdir = fullfile(prepdir, 'logs');
if ~exist(logdir, 'dir'); mkdir(logdir); end
logfile = fullfile(logdir, sprintf('%s_%s_%s_%s.log', ...
    sj_id, ses_id, 'anat', datetime('now', 'Format', 'yyyyMMdd-HHmmss')));

diary(logfile);

fprintf('\n\n\nPreprocessing: %s, %s...\n', preptype, 'anat');
fprintf('\n\n\nWorking on %s, %s...\n', sj_id, ses_id);

rawsesdir = fullfile(bidsdir, sj_id, ses_id);
rawsesanatdir = fullfile(rawsesdir, 'anat');
prepsesdir = fullfile(prepdir, sj_id, ses_id);
prepsesanatdir = fullfile(prepsesdir, 'anat');
if ~exist(prepsesanatdir, 'dir'); mkdir(prepsesanatdir); end

fsdir = fullfile(prepdir, 'freesurfer');
if ~exist(fsdir, 'dir'); mkdir(fsdir); end
cfdir = fullfile(prepdir, 'ciftify');
if ~exist(cfdir, 'dir'); mkdir(cfdir); end

metadir = fullfile(prepdir, 'metadata');
if ~exist(metadir, 'dir'); mkdir(metadir); end

%% Get a list of images

fprintf('\n\n\nGetting a list of images...\n');

anat.basename = sprintf('%s_%s', sj_id, ses_id);
anat.raw.T1 = fullfile(rawsesanatdir, '*T1w.nii.gz');
anat.raw.T1 = deblank(ls(anat.raw.T1));
anat.raw.T2 = fullfile(strrep(rawsesanatdir, ses_id, sesnext_id), '*T2w.nii.gz');
anat.raw.T2 = deblank(ls(anat.raw.T2));

%% Cross-session alignment

fprintf('\n\n\nCross-session alignment...\n');

if ses_num ~= 1
    anat.csa.T1 = fullfile(prepsesanatdir, ...
        sprintf('%s_space-anatref_T1w.nii.gz', anat.basename));
    anat.csa.ref = strrep(anat.raw.T1, ses_id, 'ses-01');
    anat.csa.xfm.os2rs = fullfile(prepsesanatdir, ...
        sprintf('%s_from-origses_to-refses.lta', anat.basename));
    anat.csa.xfmwei.os2rs = fullfile(prepsesanatdir, ...
        sprintf('%s_from-origses_to-refses.mgz', anat.basename));
    if ~exist(anat.csa.T1, 'file') || do_override
        system(sprintf('mri_robust_register --mov %s --dst %s --lta %s --mapmov %s --weights %s --iscale --satit', ...
            anat.raw.T1, anat.csa.ref, anat.csa.xfm.os2rs, anat.csa.T1, anat.csa.xfmwei.os2rs));
    end
end

%% Surface reconstruction Step 1: Bias correction & Conform (recon-all -autorecon1 -noskullstrip)

fprintf('\n\n\nSurface reconstruction Step 1: Bias correction & Conform...\n');

if ses_num == 1
    recon_all_input = anat.raw.T1;
else
    recon_all_input = anat.csa.T1;
end
if ~exist(fullfile(fsdir, anat.basename), 'dir') || do_override
    system(sprintf('recon-all -sd %s -s %s -i %s -hires -autorecon1 -noskullstrip', ...
        fsdir, anat.basename, recon_all_input));
end

%% Surface reconstruction Step 2: Skull stripping (antsBrainExtraction)

fprintf('\n\n\nSurface reconstruction Step 2: Skull stripping...\n');

brainmask = fullfile(fsdir, anat.basename, 'mri', 'brainmask.auto.mgz');
if ~exist(brainmask, 'file') || do_override
    tmpdir = tempname(tempdir); mkdir(tmpdir);
    system(sprintf('mri_convert %s %s', ...
        fullfile(fsdir, anat.basename, 'mri', 'orig.mgz'), fullfile(tmpdir, 'orig.nii.gz')));

    system(sprintf('antsBrainExtraction.sh -u 0 -d 3 -q 1 -a %s -e %s -m %s -f %s -o %s', ...
        fullfile(tmpdir, 'orig.nii.gz'), ...
        fullfile(oasisdir, 'T_template0.nii.gz'), ...
        fullfile(oasisdir, 'T_template0_BrainCerebellumProbabilityMask.nii.gz'), ...
        fullfile(oasisdir, 'T_template0_BrainCerebellumRegistrationMask.nii.gz'), ...
        fullfile(tmpdir, 'out_')));
    
    system(sprintf('mri_convert -odt uchar %s %s', ...
        fullfile(tmpdir, 'out_BrainExtractionMask.nii.gz'), fullfile(tmpdir, 'out_BrainExtractionMask.mgz')));
    system(sprintf('mri_mask %s %s %s', ...
        fullfile(fsdir, anat.basename, 'mri', 'T1.mgz'), fullfile(tmpdir, 'out_BrainExtractionMask.mgz'), brainmask));
    system(sprintf('cp %s %s', ...
        brainmask, strrep(brainmask, 'brainmask.auto.mgz', 'brainmask.mgz')));
    
    rmdir(tmpdir, 's');
end

%% Surface reconstruction Step 3: Main analysis (recon-all -autorecon2 -autorecon3)

fprintf('\n\n\nSurface reconstruction Step 3: Main analysis...\n');

wmparc = fullfile(fsdir, anat.basename, 'mri', 'wmparc.mgz');
if ~exist(wmparc, 'file') || do_override
    system(sprintf('recon-all -sd %s -s %s -T2 %s -T2pial -hires -autorecon2 -autorecon3', ...
        fsdir, anat.basename, anat.raw.T2));
end
    
%% Ciftify: ciftify_recon_all

fprintf('\n\n\nRunning ciftify_recon_all...\n');

if ~exist(fullfile(cfdir, anat.basename), 'dir') || do_override
    if ses_num == 1
        system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
            'export PATH=$PATH:%s/ciftify/bin;', ...
            'export PYTHONPATH=$PYTHONPATH:%s;', ...
            'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
            'ciftify_recon_all --ciftify-work-dir %s --fs-subjects-dir %s --resample-to-T1w32k --T2 -v %s'], ...
            cfgitdir, cfgitdir, cfgitdir, cfdir, fsdir, anat.basename));
    else
        system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
            'export PATH=$PATH:%s/ciftify/bin;', ...
            'export PYTHONPATH=$PYTHONPATH:%s;', ...
            'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
            'ciftify_recon_all --ciftify-work-dir %s --fs-subjects-dir %s --resample-to-T1w32k --read-non-lin-xfm %s --read-lin-premat %s --T2 -v %s'], ...
            cfgitdir, cfgitdir, cfgitdir, cfdir, fsdir, ...
            fullfile(cfdir, sprintf('%s_ses-01', sj_id), 'MNINonLinear', 'xfms', 'T1w2Standard_warp_noaffine.nii.gz'), ...
            fullfile(cfdir, sprintf('%s_ses-01', sj_id), 'MNINonLinear', 'xfms', 'T1w2StandardLinear.mat'), ...
            anat.basename));
    end
    system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
        'export PATH=$PATH:%s/ciftify/bin;', ...
        'export PYTHONPATH=$PYTHONPATH:%s;', ...
        'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
        'cifti_vis_recon_all snaps --ciftify-work-dir %s --verbose %s'], ...
        cfgitdir, cfgitdir, cfgitdir, cfdir, anat.basename));
end

%% Ciftify derivatives: Copy

fprintf('\n\n\nCopying ciftify derivatives...\n');

anat.xfm.fs2mni_init = fullfile(prepsesanatdir, ...
    sprintf('%s_from-fsnative_to-MNI152NLin6Asym_desc-init_xfm.mat', anat.basename));
anat.xfm.fs2mni = fullfile(prepsesanatdir, ...
    sprintf('%s_from-fsnative_to-MNI152NLin6Asym_xfm.nii.gz', anat.basename));
anat.xfm.mni2fs = fullfile(prepsesanatdir, ...
    sprintf('%s_from-MNI152NLin6Asym_to-fsnative_xfm.nii.gz', anat.basename));
anat.fs.T1prep = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-preproc_T1w.nii.gz', anat.basename));
anat.fs.T2prep = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-preproc_T2w.nii.gz', anat.basename));
anat.fs.brainmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-brain_mask.nii.gz', anat.basename));
anat.fs.T1prepmasked = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-preprocmasked_T1w.nii.gz', anat.basename));
anat.mni.T1prep = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-preproc_T1w.nii.gz', anat.basename));
anat.mni.T2prep = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-preproc_T2w.nii.gz', anat.basename));
anat.mni.brainmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-brain_mask.nii.gz', anat.basename));
anat.mni.T1prepmasked = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-preprocmasked_T1w.nii.gz', anat.basename));

if ~exist(anat.mni.T1prep, 'file') || do_override
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'xfms', 'T1w2StandardLinear.mat'), anat.xfm.fs2mni_init));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'xfms', 'T1w2Standard_warp_noaffine.nii.gz'), anat.xfm.fs2mni));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'xfms', 'Standard2T1w_warp_noaffine.nii.gz'), anat.xfm.mni2fs));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'T1w.nii.gz'), anat.fs.T1prep));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'T2w.nii.gz'), anat.fs.T2prep));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'brainmask_fs.nii.gz'), anat.fs.brainmask));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'T1w.nii.gz'), anat.mni.T1prep));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'T2w.nii.gz'), anat.mni.T2prep));
    system(sprintf('cp %s %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'brainmask_fs.nii.gz'), anat.mni.brainmask));
    system(sprintf('fslmaths %s -mas %s %s', anat.fs.T1prep, anat.fs.brainmask, anat.fs.T1prepmasked));
    system(sprintf('fslmaths %s -mas %s %s', anat.mni.T1prep, anat.mni.brainmask, anat.mni.T1prepmasked));
end

%% Ciftify derivatives: No bias-field correction (For myelin maps)

fprintf('\n\n\nGenerating T1 and T2 without bias-field corrections...\n');

anat.fs.T1orig = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_T1w.nii.gz', anat.basename));
anat.fs.T2orig = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_T2w.nii.gz', anat.basename));
anat.mni.T1orig = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_T1w.nii.gz', anat.basename));
anat.mni.T2orig = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_T2w.nii.gz', anat.basename));

if ~exist(anat.mni.T1orig, 'file') || do_override
    system(sprintf('mri_convert -odt float %s %s --conform_min -ns 1; fslreorient2std %s %s', ...
        fullfile(fsdir, anat.basename, 'mri', 'rawavg.mgz'), anat.fs.T1orig, anat.fs.T1orig, anat.fs.T1orig));
    system(sprintf('mri_convert %s %s; fslreorient2std %s %s', ...
        fullfile(fsdir, anat.basename, 'mri', 'T2.prenorm.mgz'), anat.fs.T2orig, anat.fs.T2orig, anat.fs.T2orig));
    system(sprintf('applywarp --rel --interp=spline -i %s -r %s -w %s -o %s', ...
        anat.fs.T1orig, anat.mni.T1prep, anat.xfm.fs2mni, anat.mni.T1orig));
    system(sprintf('fslmaths %s -max 0 %s', anat.mni.T1orig, anat.mni.T1orig));
    system(sprintf('applywarp --rel --interp=spline -i %s -r %s -w %s -o %s', ...
        anat.fs.T2orig, anat.mni.T1prep, anat.xfm.fs2mni, anat.mni.T2orig));
    system(sprintf('fslmaths %s -max 0 %s', anat.mni.T2orig, anat.mni.T2orig));
end

%% Ciftify derivatives: GM/WM/CSF masks

fprintf('\n\n\nGenerating GM/WM/CSF masks...\n');

anat.fs.GMmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-GM_mask.nii.gz', anat.basename));
anat.fs.WMmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-WM_mask.nii.gz', anat.basename));
anat.fs.CSFmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-CSF_mask.nii.gz', anat.basename));
anat.mni.GMmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-GM_mask.nii.gz', anat.basename));
anat.mni.WMmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-WM_mask.nii.gz', anat.basename));
anat.mni.CSFmask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-CSF_mask.nii.gz', anat.basename));

if ~exist(anat.mni.GMmask, 'file') || do_override
    system(sprintf('mri_binarize --i %s --gm --o %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'aparc+aseg.nii.gz'), anat.fs.GMmask));
    system(sprintf('mri_binarize --i %s --all-wm --o %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'aparc+aseg.nii.gz'), anat.fs.WMmask));
    system(sprintf('mri_binarize --i %s --ventricles --match 15, 24 --o %s', ...
        fullfile(cfdir, anat.basename, 'T1w', 'aparc+aseg.nii.gz'), anat.fs.CSFmask));
    system(sprintf('mri_binarize --i %s --gm --o %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'aparc+aseg.nii.gz'), anat.mni.GMmask));
    system(sprintf('mri_binarize --i %s --all-wm --o %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'aparc+aseg.nii.gz'), anat.mni.WMmask));
    system(sprintf('mri_binarize --i %s --ventricles --match 15, 24 --o %s', ...
        fullfile(cfdir, anat.basename, 'MNINonLinear', 'aparc+aseg.nii.gz'), anat.mni.CSFmask));
end

%% Ciftify derivatives: Erode WM/CSF masks

fprintf('\n\n\nEroding WM/CSF masks...\n');

anat.fs.WMeromask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-WMero4_mask.nii.gz', anat.basename));
anat.fs.CSFeromask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-fsnative_desc-CSFero2_mask.nii.gz', anat.basename));
anat.mni.WMeromask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-WMero4_mask.nii.gz', anat.basename));
anat.mni.CSFeromask = fullfile(prepsesanatdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-CSFero2_mask.nii.gz', anat.basename));

if ~exist(anat.mni.WMeromask, 'file') || do_override
    system(sprintf('fslmaths %s -ero -ero -ero -ero %s', anat.fs.WMmask, anat.fs.WMeromask));
    system(sprintf('fslmaths %s -ero -ero %s', anat.fs.CSFmask, anat.fs.CSFeromask));
    system(sprintf('applywarp --rel --interp=nn -i %s -r %s -w %s -o %s', ...
        anat.fs.WMeromask, anat.mni.T1prep, anat.xfm.fs2mni, anat.mni.WMeromask));
    system(sprintf('applywarp --rel --interp=nn -i %s -r %s -w %s -o %s', ...
        anat.fs.CSFeromask, anat.mni.T1prep, anat.xfm.fs2mni, anat.mni.CSFeromask));
end

%% Save metadata and finish

fprintf('\n\n\nSaving metadata...\n');

metafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, ses_id, 'anat'));
if ~exist(metafile, 'file') || do_override
    save(metafile, 'anat');
end

fprintf('\n\n\nDone.\n');

diary off;

end
