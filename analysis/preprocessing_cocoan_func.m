function preprocessing_cocoan_func(sj_num, ses_num, run_num, cfgitdir, varargin)

setenv('DYLD_LIBRARY_PATH', [getenv('FREESURFER_HOME') '/lib/gcc/lib' ':/opt/X11/lib/flat_namespace']);

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

tasklabels = {'task-rest', 'task-rating_run-1', ...
    'task-rating_run-2', 'task-rating_run-3', ...
    'task-speak_run-1', 'task-speak_run-2', ...
    'task-listen_run-1', 'task-listen_run-2'};
tasklabel = tasklabels{run_num};

logdir = fullfile(prepdir, 'logs');
if ~exist(logdir, 'dir'); mkdir(logdir); end
logfile = fullfile(logdir, sprintf('%s_%s_%s_%s_%s.log', ...
    sj_id, ses_id, 'func', tasklabel, datetime('now', 'Format', 'yyyyMMdd-HHmmss')));

diary(logfile);

fprintf('\n\n\nPreprocessing: %s, %s...\n', preptype, 'func');
fprintf('\n\n\nWorking on %s, %s, %s...\n', sj_id, ses_id, tasklabel);

rawsesdir = fullfile(bidsdir, sj_id, ses_id);
rawsesfuncdir = fullfile(rawsesdir, 'func');
rawsesfmapdir = fullfile(rawsesdir, 'fmap');
prepsesdir = fullfile(prepdir, sj_id, ses_id);
prepsesfuncdir = fullfile(prepsesdir, 'func');
if ~exist(prepsesfuncdir, 'dir'); mkdir(prepsesfuncdir); end
prepsesfmapdir = fullfile(prepsesdir, 'fmap');
if ~exist(prepsesfmapdir, 'dir'); mkdir(prepsesfmapdir); end

cfdir = fullfile(prepdir, 'ciftify');

metadir = fullfile(prepdir, 'metadata');
anatmetafile = fullfile(metadir, sprintf('%s_%s_%s_metadata.mat', sj_id, 'ses-01', 'anat'));
load(anatmetafile, 'anat');

%% Get a list of images

fprintf('\n\n\nGetting a list of images...\n');

func.basename = sprintf('%s_%s_%s', sj_id, ses_id, tasklabel);
func.raw.bold = fullfile(rawsesfuncdir, sprintf('%s*_bold.nii.gz', func.basename));
func.raw.bold = deblank(ls(func.raw.bold));
func.raw.json = strrep(func.raw.bold, 'nii.gz', 'json');
func.raw.jsondat = jsondecode(fileread(func.raw.json));
fmap.basename = sprintf('%s_%s_acq-%s', sj_id, ses_id, 'func');
fmap.raw.epi = fullfile(rawsesfmapdir, sprintf('%s*_epi.nii.gz', fmap.basename));
fmap.raw.epi = split(deblank(ls(fmap.raw.epi)));
fmap.raw.json = strrep(fmap.raw.epi, 'nii.gz', 'json');
for i = 1:numel(fmap.raw.json)
    fmap.raw.jsondat{i,1} = jsondecode(fileread(fmap.raw.json{i}));
end

%% Discard first 10 secs

fprintf('\n\n\nRemove initial dummy volumes (first 10 secs)...\n');

func.rmd.bold = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-removedummy_bold.nii.gz', func.basename));

func.rmd.tdummy = 10; % 10 seconds
func.rmd.ndummy = ceil(func.rmd.tdummy ./ func.raw.jsondat.RepetitionTime); % 22 images for TR 0.46 s
[~, sysout] = system(sprintf('fslnvols %s', func.raw.bold));
func.rmd.nvol = str2double(sysout) - func.rmd.ndummy;

if ~exist(func.rmd.bold, 'file') || do_override
    system(sprintf('fslroi %s %s %d %d', func.raw.bold, func.rmd.bold, func.rmd.ndummy, func.rmd.nvol));
    system(sprintf('fslmaths %s %s -odt float', func.rmd.bold, func.rmd.bold));
end

%% Motion correction

fprintf('\n\n\nMotion correction...\n');

func.mc.bold = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-mc_bold.nii.gz', func.basename));
func.mc.rp = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-mc_rp.txt', func.basename));
func.mc.boldmean = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-mcmean_bold.nii.gz', func.basename));

if ~exist(func.mc.bold, 'file') || do_override
    system(sprintf('3dvolreg -Fourier -twopass -prefix %s -1Dfile %s %s', ...
        func.mc.bold, func.mc.rp, func.rmd.bold));
    system(sprintf('fslmaths %s -Tmean %s', func.mc.bold, func.mc.boldmean));
end

%% Merge fieldmap and get acquisition parameters

fprintf('\n\n\nMerge fieldmap images...\n');

fmap.merged.epi = fullfile(prepsesfmapdir, ...
    sprintf('%s_desc-merged_epi.nii.gz', fmap.basename));
if ~exist(fmap.merged.epi, 'file') || do_override
    system(sprintf(['fslmerge -t %s' repmat(' %s', 1, numel(fmap.raw.epi))], ...
        fmap.merged.epi, fmap.raw.epi{:}));
    system(sprintf('fslmaths %s %s -odt float', fmap.merged.epi, fmap.merged.epi));
end

fmap.merged.acqparam = fullfile(prepsesfmapdir, ...
    sprintf('%s_desc-merged_acqparam.txt', fmap.basename));
fmap.merged.acqparamdat = [];
fmap.merged.volinfo = [];
for i = 1:numel(fmap.raw.epi)
    pedir = [1 -1] * strcmp(fmap.raw.jsondat{i}.PhaseEncodingDirection, {'i', 'j', 'k'; 'i-', 'j-', 'k-'});
    [~, sysout] = system(sprintf('fslinfo %s | grep dim | grep -v pix | awk ''{print $2}''', fmap.raw.epi{i}));
    fmap.merged.volinfo(i,:) = str2num(sysout);
    fmap.merged.acqparamdat = [fmap.merged.acqparamdat; ...
        repmat([pedir, fmap.raw.jsondat{i}.TotalReadoutTime], fmap.merged.volinfo(i,4), 1)];
end
if ~exist(fmap.merged.acqparam, 'file') || do_override
    writematrix(fmap.merged.acqparamdat, fmap.merged.acqparam, 'Delimiter', '\t');
end

%% Coregister fieldmap images

fprintf('\n\n\nCoregister fieldmap images to functional images...\n');

fmap.coreg.funcidx = [];
for i = 1:numel(fmap.raw.epi)
    if strcmp(func.raw.jsondat.PhaseEncodingDirection, fmap.raw.jsondat{i}.PhaseEncodingDirection)
        fmap.coreg.funcidx = [fmap.coreg.funcidx, i];
    end
end
fmap.coreg.xfm.fmap2func = fullfile(prepsesfmapdir, ...
    sprintf('%s_from-fmap_to-func_%s_xfm.mat', fmap.basename, tasklabel));
if ~exist(fmap.coreg.xfm.fmap2func, 'file') || do_override
    system(sprintf('flirt -in %s -ref %s -omat %s -cost mutualinfo -dof 6', ...
        fmap.raw.epi{fmap.coreg.funcidx(1)}, func.mc.boldmean, fmap.coreg.xfm.fmap2func)); % internally use 1st vol
end

fmap.coreg.epi = fullfile(prepsesfmapdir, ...
    sprintf('%s_space-func_desc-merged_%s_epi.nii.gz', fmap.basename, tasklabel));
if ~exist(fmap.coreg.epi, 'file') || do_override
    system(sprintf('flirt -in %s -ref %s -applyxfm -init %s -interp spline -out %s', ...
        fmap.merged.epi, func.mc.boldmean, fmap.coreg.xfm.fmap2func, fmap.coreg.epi));
    system(sprintf('fslmaths %s -max 0 %s', fmap.coreg.epi, fmap.coreg.epi));
end

%% Topup

fprintf('\n\n\nTopup: Estimating susceptibility distortion...\n');

fmap.topup.basename = fullfile(prepsesfmapdir, ...
    sprintf('%s_space-func_desc-preproc_%s', fmap.basename, tasklabel));
fmap.topup.epi = fullfile(prepsesfmapdir, ...
    sprintf('%s_space-func_desc-preproc_%s_epi.nii.gz', fmap.basename, tasklabel));

if any(rem(fmap.merged.volinfo(:,1:3), 2) == 1, 1:2)
    fmap.topup.config = fullfile(getenv('FSLDIR'), 'etc/flirtsch/b02b0_1.cnf'); % subsampling 1, takes longer
else
    fmap.topup.config = fullfile(getenv('FSLDIR'), 'etc/flirtsch/b02b0.cnf');
end

if ~exist(fmap.topup.epi, 'file') || do_override
    system(sprintf('topup --imain=%s --datain=%s --config=%s --out=%s --iout=%s --verbose', ...
        fmap.coreg.epi, fmap.merged.acqparam, fmap.topup.config, fmap.topup.basename, fmap.topup.epi));
end

%% Applytopup

fprintf('\n\n\nApplytopup: Correcting susceptibility distortion...\n');

func.dc.bold = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-dc_bold.nii.gz', func.basename));
func.dc.boldmean = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-dcmean_bold.nii.gz', func.basename));
[~, func.dc.inindex] = find(repelem(1:size(fmap.merged.volinfo,1), 1, fmap.merged.volinfo(:,4)) == fmap.coreg.funcidx(1));

if ~exist(func.dc.bold, 'file') || do_override
    system(sprintf('applytopup --imain=%s --inindex=%d --topup=%s --datain=%s --method=jac --out=%s --verbose', ...
        func.mc.bold, func.dc.inindex(1), fmap.topup.basename, fmap.merged.acqparam, func.dc.bold));
    system(sprintf('fslmaths %s -max 0 %s', func.dc.bold, func.dc.bold));
    system(sprintf('fslmaths %s -Tmean %s', func.dc.bold, func.dc.boldmean));
end

%% Cross-run alignment

fprintf('\n\n\nCross-run alignment...\n');

if run_num ~= 1
    func.cra.ref = strrep(func.dc.boldmean, tasklabel, tasklabels{1});
    func.cra.xfm.or2rr = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-origrun_to-refrun_xfm.mat', func.basename));
    func.cra.xfm.rr2or = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-refrun_to-origrun_xfm.mat', func.basename));
    if ~exist(func.cra.xfm.or2rr, 'file') || do_override
        system(sprintf('flirt -in %s -ref %s -omat %s -cost normcorr -dof 6', ...
            func.dc.boldmean, func.cra.ref, func.cra.xfm.or2rr));
        system(sprintf('convert_xfm -omat %s -inverse %s', func.cra.xfm.rr2or, func.cra.xfm.or2rr));
    end
end

%% Cross-session alignment

fprintf('\n\n\nCross-session alignment...\n');

if ses_num ~= 1 && run_num == 1
    func.csa.ref = strrep(func.dc.boldmean, ses_id, 'ses-01');
    func.csa.xfm.os2rs = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-origses_to-refses_xfm.mat', func.basename));
    func.csa.xfm.rs2os = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-refses_to-origses_xfm.mat', func.basename));
    if ~exist(func.csa.xfm.os2rs, 'file') || do_override
        system(sprintf('flirt -in %s -ref %s -omat %s -cost normcorr -dof 6', ...
            func.dc.boldmean, func.csa.ref, func.csa.xfm.os2rs));
        system(sprintf('convert_xfm -omat %s -inverse %s', func.csa.xfm.rs2os, func.csa.xfm.os2rs));
    end
end

%% Coregistration

fprintf('\n\n\nCoregister functional images to structural images...\n');

func.coreg.xfm.func2fs = fullfile(prepsesfuncdir, ...
    sprintf('%s_from-func_to-fsnative_xfm.mat', func.basename));
func.coreg.xfm.fs2func = fullfile(prepsesfuncdir, ...
    sprintf('%s_from-fsnative_to-func_xfm.mat', func.basename));

if ses_num == 1 && run_num == 1
    
    func.coreg.xfm.func2fs_initfslflirt = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-func_to-fsnative_desc-initfslflirt_xfm.mat', func.basename));
    func.coreg.xfm.func2fs_initfslbbr = fullfile(prepsesfuncdir, ...
        sprintf('%s_from-func_to-fsnative_desc-initfslbbr_xfm.mat', func.basename));
    if ~exist(func.coreg.xfm.func2fs, 'file') || do_override
        system(sprintf('flirt -in %s -ref %s -omat %s -cost mutualinfo -dof 6', ...
            func.dc.boldmean, anat.fs.T1prepmasked, func.coreg.xfm.func2fs_initfslflirt));
        system(sprintf('flirt -in %s -ref %s -init %s -omat %s -cost bbr -dof 6 -schedule %s -wmseg %s', ...
            func.dc.boldmean, anat.fs.T1prep, func.coreg.xfm.func2fs_initfslflirt, func.coreg.xfm.func2fs_initfslbbr, ...
            fullfile(getenv('FSLDIR'), 'etc/flirtsch/bbr.sch'), anat.fs.WMmask));
        
        tmpdir = tempname(tempdir); mkdir(tmpdir);
        t1mgz = fullfile(prepdir, 'freesurfer', anat.basename, 'mri', 'T1.mgz');
        system(sprintf('mri_convert %s %s', t1mgz, fullfile(tmpdir, 'T1.nii.gz')));
        system(sprintf('fslreorient2std -m %s %s', fullfile(tmpdir, 'fs2fslfs.mat'), fullfile(tmpdir, 'T1.nii.gz')));
        system(sprintf('convert_xfm -omat %s -inverse %s', fullfile(tmpdir, 'fslfs2fs.mat'), fullfile(tmpdir, 'fs2fslfs.mat')));
        system(sprintf('convert_xfm -omat %s -concat %s %s', fullfile(tmpdir, 'bbr.mat'), fullfile(tmpdir, 'fslfs2fs.mat'), func.coreg.xfm.func2fs_initfslbbr));
        system(sprintf('lta_convert --infsl %s --outlta %s --src %s --trg %s --subject %s', ...
            fullfile(tmpdir, 'bbr.mat'), fullfile(tmpdir, 'bbr.lta'), func.dc.boldmean, t1mgz, anat.basename));

        system(sprintf('export SUBJECTS_DIR=%s; bbregister --s %s --mov %s --bold --init-reg %s --reg %s --fslmat %s', ...
            fullfile(prepdir, 'freesurfer'), anat.basename, func.dc.boldmean, ...
            fullfile(tmpdir, 'bbr.lta'), fullfile(tmpdir, 'bbreg.lta'), fullfile(tmpdir, 'bbreg.mat')));
        system(sprintf('convert_xfm -omat %s -concat %s %s', func.coreg.xfm.func2fs, fullfile(tmpdir, 'fs2fslfs.mat'), fullfile(tmpdir, 'bbreg.mat')));
        system(sprintf('convert_xfm -omat %s -inverse %s', func.coreg.xfm.fs2func, func.coreg.xfm.func2fs));
        
        rmdir(tmpdir, 's');
    end
    
elseif ses_num ~= 1 && run_num == 1
    
    if ~exist(func.coreg.xfm.func2fs, 'file') || do_override
        system(sprintf('convert_xfm -omat %s -concat %s %s', ...
            func.coreg.xfm.func2fs, ...
            strrep(func.coreg.xfm.func2fs, ses_id, 'ses-01'), ...
            func.csa.xfm.os2rs));
        system(sprintf('convert_xfm -omat %s -inverse %s', func.coreg.xfm.fs2func, func.coreg.xfm.func2fs));
    end
 
elseif run_num ~= 1
    
    if ~exist(func.coreg.xfm.func2fs, 'file') || do_override
        system(sprintf('convert_xfm -omat %s -concat %s %s', ...
            func.coreg.xfm.func2fs, ...
            strrep(func.coreg.xfm.func2fs, tasklabel, tasklabels{1}), ...
            func.cra.xfm.or2rr));
        system(sprintf('convert_xfm -omat %s -inverse %s', func.coreg.xfm.fs2func, func.coreg.xfm.func2fs));
    end
       
end

%% Normalization

fprintf('\n\n\nNormalization to MNI template...\n');

func.mni.bold = fullfile(prepsesfuncdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz', func.basename));
func.mni.boldmean = fullfile(prepsesfuncdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-preprocmean_bold.nii.gz', func.basename));

if ~exist(func.mni.bold, 'file') || do_override
    tmpdir = tempname(tempdir); mkdir(tmpdir);
    system(sprintf('fslsplit %s %s', func.dc.bold, fullfile(tmpdir, 'pre')));
    system(sprintf('for i in {0..%d}; do applywarp --rel --interp=spline -i %s -r %s --premat=%s -w %s -o %s; done', ...
        func.rmd.nvol-1, fullfile(tmpdir, 'pre$(printf "%.4d" $i).nii.gz'), anat.mni.T1prep, ...
        func.coreg.xfm.func2fs, anat.xfm.fs2mni, fullfile(tmpdir, 'post$(printf "%.4d" $i).nii.gz')));
    system(sprintf('cd %s; fslmerge -t %s %s', tmpdir, func.mni.bold, 'post*'));
    system(sprintf('fslmaths %s -Tmean %s', func.mni.bold, func.mni.boldmean));
    rmdir(tmpdir, 's');
end

%% Denoising (Linear trend, 6RP, aCompCor 5WM/5CSF, Censoring with LPF(0.1Hz)-FD > 0.15mm, BPF 0.005-0.1Hz OR HPF 0.005Hz)

fprintf('\n\n\nDenoising...\n');

func.den.bold = fullfile(prepsesfuncdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-denoise_bold.nii.gz', func.basename));
func.den.boldmean = fullfile(prepsesfuncdir, ...
    sprintf('%s_space-MNI152NLin6Asym_res-2_desc-denoisemean_bold.nii.gz', func.basename));
func.den.nui = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-confounds_timeseries.txt', func.basename));
func.den.cen = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-censoring_timeseries.txt', func.basename));
func.den.cenlpfthr = 0.1;
func.den.cenfdthr = 0.15;
func.den.npoly = 1;
if strcmp(tasklabel, 'task-rest')
    func.den.bpf = [0.005 0.1];
else
    func.den.bpf = [0.005 9999];
end
func.den.ort = fullfile(prepsesfuncdir, ...
    sprintf('%s_desc-ort_timeseries.txt', func.basename));

if ~exist(func.den.bold, 'file') || do_override
    rp = importdata(func.mc.rp);
    
    mnibolddat = read_avw(func.mni.bold);
    mnibolddat = reshape(mnibolddat, [], size(mnibolddat,4));
    WMeromaskdat = read_avw(anat.mni.WMeromask);
    CSFeromaskdat = read_avw(anat.mni.CSFeromask);
    WMdat = mnibolddat(logical(WMeromaskdat(:)), :);
    CSFdat = mnibolddat(logical(CSFeromaskdat(:)), :);
    mnibolddat = [];
    [~, WMpc] = pca(detrend(WMdat.'));
    [~, CSFpc] = pca(detrend(CSFdat.'));
    
    writematrix([rp WMpc(:,1:5) CSFpc(:,1:5)], func.den.nui, 'Delimiter', '\t');
    
    [B, A] = butter(2, func.den.cenlpfthr / ((1/func.raw.jsondat.RepetitionTime) / 2), 'low');
    filtrp = filtfilt(B, A, rp);
    filtrp = [filtrp(:,4:6), filtrp(:,1:3) * pi / 180 * 50];
    filtfd = sum(abs(diff(filtrp([1 1:end], :))), 2);
    writematrix(filtfd <= func.den.cenfdthr, func.den.cen);
    
    tmpdir = tempname(tempdir); mkdir(tmpdir);
    system(sprintf('3dTproject -input %s -censor %s -cenmode NTRP -polort %d -ort %s -passband %f %f -TR %f -prefix %s -verb', ...
        func.mni.bold, func.den.cen, func.den.npoly, func.den.nui, func.den.bpf, func.raw.jsondat.RepetitionTime, fullfile(tmpdir, 'den.nii.gz')));
    system(sprintf('mv %s %s', fullfile(tmpdir, 'den.nii.gz'), func.den.bold));
    system(sprintf('mv %s %s', fullfile(tmpdir, 'den.nii.gz.ort.1D'), func.den.ort));
    system(sprintf('fslmaths %s -add %s %s', func.den.bold, func.mni.boldmean, func.den.bold));
    system(sprintf('fslmaths %s -Tmean %s', func.den.bold, func.den.boldmean));
    
    rmdir(tmpdir, 's');
end

%% Ciftify: ciftify_subject_fmri

fprintf('\n\n\nRunning ciftify_subject_fmri...\n');

[~, func.cf.mni.basename] = fileparts(strrep(func.mni.bold, '.nii.gz', ''));
func.cf.mni.fwhm = 6;
func.cf.mni.rbold = fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.mni.basename, ...
    sprintf('%s_Atlas_s%d.dtseries.nii', func.cf.mni.basename, 0));
func.cf.mni.sbold = fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.mni.basename, ...
    sprintf('%s_Atlas_s%d.dtseries.nii', func.cf.mni.basename, func.cf.mni.fwhm));
if ~exist(func.cf.mni.rbold, 'file') || do_override
    system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
        'export PATH=$PATH:%s/ciftify/bin;', ...
        'export PYTHONPATH=$PYTHONPATH:%s;', ...
        'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
        'ciftify_subject_fmri --ciftify-work-dir %s --SmoothingFWHM %d ', ...
        '--func-ref %s --already-in-MNI --OutputSurfDiagnostics -v %s %s %s'], ...
        cfgitdir, cfgitdir, cfgitdir, cfdir, func.cf.mni.fwhm, func.mni.boldmean, func.mni.bold, anat.basename, func.cf.mni.basename));
    system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
        'export PATH=$PATH:%s/ciftify/bin;', ...
        'export PYTHONPATH=$PYTHONPATH:%s;', ...
        'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
        'cifti_vis_fmri snaps --ciftify-work-dir %s --SmoothingFWHM %d --meanfunc %s --verbose %s %s'], ...
        cfgitdir, cfgitdir, cfgitdir, cfdir, func.cf.mni.fwhm, func.mni.boldmean, func.cf.mni.basename, anat.basename));
    system(sprintf('rm %s.nii.gz', fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.mni.basename, func.cf.mni.basename)));
end

[~, func.cf.den.basename] = fileparts(strrep(func.den.bold, '.nii.gz', ''));
func.cf.den.fwhm = 6;
func.cf.den.rbold = fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.den.basename, ...
    sprintf('%s_Atlas_s%d.dtseries.nii', func.cf.den.basename, 0));
func.cf.den.sbold = fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.den.basename, ...
    sprintf('%s_Atlas_s%d.dtseries.nii', func.cf.den.basename, func.cf.den.fwhm));
if ~exist(func.cf.den.rbold, 'file') || do_override
    system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
        'export PATH=$PATH:%s/ciftify/bin;', ...
        'export PYTHONPATH=$PYTHONPATH:%s;', ...
        'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
        'ciftify_subject_fmri --ciftify-work-dir %s --SmoothingFWHM %d ', ...
        '--func-ref %s --already-in-MNI --OutputSurfDiagnostics -v %s %s %s'], ...
        cfgitdir, cfgitdir, cfgitdir, cfdir, func.cf.den.fwhm, func.den.boldmean, func.den.bold, anat.basename, func.cf.den.basename));
    system(sprintf(['source ~/.ciftifyvenv/bin/activate;', ...
        'export PATH=$PATH:%s/ciftify/bin;', ...
        'export PYTHONPATH=$PYTHONPATH:%s;', ...
        'export CIFTIFY_TEMPLATES=%s/ciftify/data;', ...
        'cifti_vis_fmri snaps --ciftify-work-dir %s --SmoothingFWHM %d --meanfunc %s --verbose %s %s'], ...
        cfgitdir, cfgitdir, cfgitdir, cfdir, func.cf.den.fwhm, func.den.boldmean, func.cf.den.basename, anat.basename));
    system(sprintf('rm %s.nii.gz', fullfile(cfdir, anat.basename, 'MNINonLinear', 'Results', func.cf.den.basename, func.cf.den.basename)));
end

%% Save metadata and finish

metafile = fullfile(metadir, sprintf('%s_%s_%s_%s_metadata.mat', sj_id, ses_id, 'func', tasklabel));
if ~exist(metafile, 'file') || do_override
    save(metafile, 'func', 'fmap');
end

fprintf('\n\n\nDone.\n');
diary off;

end