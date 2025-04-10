%% 

maxNumCompThreads(56);
bidsdir = '/Volumes/cocoanlab02/data/DEIPP/imaging/bids_dataset';
rscdir = '/Volumes/cocoanlab01/resources';
gitdir = '/Users/jaejoong/github';
sj_num = 1;
% sj_num = 2;
% sj_num = 3;
opt = struct;
opt.dfc = 'ets';
opt.parc = 'indparc';
% opt.parc = 'schaefer';
opt.bintype = 'rbin';
opt.nbin = 10;
opt.ntrain = 23;
% opt.ntrain = 28;
% opt.ntrain = 13;
opt.hout = 'none';
opt.algorithm = 'LASSOOLSPCR';
opt.param = fliplr(1:3*opt.nbin*opt.ntrain-2);

%%

analysis_resting_func_net_indparc(sj_num, fullfile(rscdir, 'fmri_toolboxes/BCT/2019_03_03_BCT'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    'bidsdir', bidsdir);

%%

analysis_resting_func_net_schaefer(sj_num, fullfile(rscdir, 'fmri_toolboxes/BCT/2019_03_03_BCT'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    fullfile(gitdir, 'Schaefer2018_LocalGlobal_Parcellations_HCP/fslr32k/cifti/Schaefer2018_500Parcels_7Networks_order.dlabel.nii'), ...
    fullfile(gitdir, 'Schaefer2018_LocalGlobal_Parcellations_HCP/fslr32k/cifti/Schaefer2018_500Parcels_7Networks_order_info.txt'), ...
    'bidsdir', bidsdir);

%%

analysis_predict_rating_func_observation(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_func_headmov(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_func_feature(sj_num, opt, 'bidsdir', bidsdir);

%%

analysis_predict_rating_func_modeling(sj_num, opt, 'bidsdir', bidsdir);

%%

analysis_predict_rating_func_test(sj_num, opt, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_results(sj_num, opt, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_headmovresults(sj_num, opt, 'bidsdir', bidsdir);

%% 

analysis_predict_rating_func_trainsize(sj_num, opt, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_trainsize(sj_num, opt, 'bidsdir', bidsdir);

%%

analysis_predict_rating_func_featimp(sj_num, opt, 'bidsdir', bidsdir);

%%

analysis_visualize_model_weights(sj_num, opt, fullfile(rscdir, 'spm12'), 'bidsdir', bidsdir);

%%

maxNumCompThreads(56);
bidsdir = '/Volumes/cocoanlab02/data/DEIPP/imaging/bids_dataset';

tesj_num = 1;
teopt = struct;
teopt.dfc = 'ets';
teopt.parc = 'indparc';
% teopt.parc = 'schaefer';
teopt.bintype = 'rbin';
teopt.nbin = 10;
teopt.ntrain = 23;
teopt.hout = 'none';
teopt.algorithm = 'LASSOOLSPCR';
teopt.param = fliplr(1:3*teopt.nbin*teopt.ntrain-2);

trsj_num = 2;
tropt = struct;
tropt.dfc = 'ets';
tropt.parc = 'indparc';
% tropt.parc = 'schaefer';
tropt.bintype = 'rbin';
tropt.nbin = 10;
tropt.ntrain = 28;
tropt.hout = 'none';
tropt.algorithm = 'LASSOOLSPCR';
tropt.param = fliplr(1:3*tropt.nbin*tropt.ntrain-2);

%%

analysis_predict_rating_func_crosstest(tesj_num, teopt, trsj_num, tropt, 'bidsdir', bidsdir);
analysis_predict_rating_func_crosstest(trsj_num, tropt, tesj_num, teopt, 'bidsdir', bidsdir);

%%

analysis_visualize_model_crossresults(tesj_num, trsj_num, tropt, 'bidsdir', bidsdir);
analysis_visualize_model_crossresults(trsj_num, tesj_num, teopt, 'bidsdir', bidsdir);
