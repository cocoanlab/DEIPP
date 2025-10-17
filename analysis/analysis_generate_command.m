%% 

maxNumCompThreads(56);
bidsdir = '/Volumes/cocoanlab02/data/DEIPP/imaging/bids_dataset';
rscdir = '/Volumes/cocoanlab01/resources';
gitdir = '/Users/jaejoong/github';
mricrogldir = '/Applications/MRIcroGL.app/Contents/MacOS/MRIcroGL';
sj_num = 1;
% sj_num = 2;
% sj_num = 3;

%%

analysis_MSC(sj_num, fullfile(gitdir, 'MSCcodebase'), fullfile(gitdir, 'Infomap'), 'bidsdir', bidsdir);

%%

analysis_SCAN(sj_num, fullfile(gitdir, 'MSCcodebase'), fullfile(gitdir, 'Infomap'), 'bidsdir', bidsdir);

%%

analysis_parcellation(sj_num, ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    'bidsdir', bidsdir);

%%

analysis_visualize_parcellation(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%%

analysis_rating(sj_num, 'bidsdir', bidsdir);

%%

analysis_visualize_rating(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_feature(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_modeling(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_test(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_test(sj_num, 'bidsdir', bidsdir);

%%

analysis_parcellation(sj_num, ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    'bidsdir', bidsdir, 'parctype', 'schaefer', ...
    'ctxparcname', fullfile(gitdir, 'Schaefer2018_LocalGlobal_Parcellations_HCP/fslr32k/cifti/Schaefer2018_500Parcels_7Networks_order.dlabel.nii'), ...
    'ctxparclabelname', fullfile(gitdir, 'Schaefer2018_LocalGlobal_Parcellations_HCP/fslr32k/cifti/Schaefer2018_500Parcels_7Networks_order_info.txt'));

%%

analysis_predict_rating_feature(sj_num, 'bidsdir', bidsdir, 'parctype', 'schaefer');

%%

analysis_predict_rating_modeling(sj_num, 'bidsdir', bidsdir, 'parctype', 'schaefer');

%%

analysis_predict_rating_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'schaefer');

%% 

analysis_visualize_model_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'schaefer');

%%

analysis_parcellation(sj_num, ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    'bidsdir', bidsdir, 'parctype', 'indparcwvis');

%%

analysis_predict_rating_feature(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparcwvis');

%%

analysis_predict_rating_modeling(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparcwvis');

%%

analysis_predict_rating_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparcwvis');

%% 

analysis_visualize_model_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparcwvis');

%%

analysis_parcellation(sj_num, ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'), ...
    fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'), ...
    'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%%

analysis_predict_rating_feature(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%%

analysis_predict_rating_modeling(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%%

analysis_predict_rating_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%% 

analysis_visualize_model_test(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%%

analysis_predict_rating_headmov(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_headmov(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_testonrun1(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_testonrest(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_testonrest(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_testonrest(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%% 

analysis_visualize_model_testonrest(sj_num, 'bidsdir', bidsdir, 'parctype', 'indparconlyvis');

%%

analysis_predict_rating_ToPStest(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_ToPS(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_modelingses12(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_ses12test(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_ses12test(sj_num, 'bidsdir', bidsdir);

%% 

analysis_predict_rating_trainsize(sj_num, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_trainsize(sj_num, 'bidsdir', bidsdir);

%%

analysis_predict_rating_featimp(sj_num, 'bidsdir', bidsdir);

%%

analysis_visualize_model_featimp(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_weights(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%% 

analysis_visualize_model_netimpwei(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%%

analysis_univariate_rating_firstlevel(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%%

analysis_univariate_rating_secondlevel(sj_num, fullfile(rscdir, 'spm12'), mricrogldir, 'bidsdir', bidsdir);

%%

analysis_predict_rating_crosstest(1, 2, 'bidsdir', bidsdir);

%%

analysis_predict_rating_crosstest(2, 1, 'bidsdir', bidsdir);

%%

analysis_visualize_model_crosstest(1, 2, 'bidsdir', bidsdir);

%%

analysis_visualize_model_crosstest(2, 1, 'bidsdir', bidsdir);

