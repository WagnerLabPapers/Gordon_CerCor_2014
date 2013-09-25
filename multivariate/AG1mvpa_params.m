function S= AG_mvpa_params(subj_id)
%
% generate classification parameters
% Alan Gordon, Stanford University

par = AG1Params(subj_id);
S.par = par;
S.exp_name = 'AG';

%% File Location Parameters

S.subj_id = subj_id;
S.expt_dir = '/Users/gordonam/Studies/AG1/fmri_data2/';
S.mvpa_dir = [S.expt_dir  S.subj_id '/mvpa'];
S.anat_dir = [S.expt_dir  S.subj_id '/anat'];
S.group_mvpa_dir = [S.expt_dir 'mvpa_files'];
S.univar_dir = [S.expt_dir S.subj_id '/' 'EncAnalysis'];
S.workspace_dir = [S.expt_dir S.subj_id '/' 'mvpa'];

%% Condition Parameters

S.trainTask = 'enc';
S.testTask = 'enc';

if strcmp(S.trainTask, S.testTask)
    %train and test on same sets of data
    S.sameTaskTrainAndTest = 1;
    S.xval = 1;
else
    %train and test on independent sets of data
    S.sameTaskTrainAndTest = 0;
    S.xval = 0;
end

if S.trainTask == 'enc'
    S.onsetsTrainDir = [S.expt_dir subj_id '/EncAnalysis/'];
    S.condsTrain = {{'PersonConfCor'}  {'SceneConfCor'}} ;
    S.TrainRuns = [1:6];
elseif S.trainTask == 'ret'
    S.onsetsTrainDir = [S.expt_dir subj_id '/RetAnalysis/'];
    S.condsTrain = {{'PCorConf'} {'SCorConf'}};
    S.TrainRuns = [7:12];
end
    
if S.testTask == 'enc'
    S.onsetsTestDir = [S.expt_dir subj_id '/EncAnalysis/'];
    S.condsTest = {{'PersonConfCor'}  {'SceneConfCor'}} ;
    S.TestRuns = [1:6];
elseif S.testTask == 'ret'
    S.onsetsTestDir = [S.expt_dir subj_id '/RetAnalysis/'];
    S.condsTest = {{'PIncConf'} {'SIncConf'}};
    S.TestRuns = [7:12];
end

% duration
S.durTrain = sum(par.numvols(S.TrainRuns)) * par.TR;
S.durTest = sum(par.numvols(S.TestRuns)) * par.TR;

% types of importance maps
S.impType = {'pos' 'neg' 'both' 'raw'};

S.nFolds = 10; %1 for leave-one-out xvalidation

% what cross-validation scheme do we use?
if S.xval
    S.thisSelector = 'randomNFold_xval';
    S.thisSigIntenseSelector = 'randomNFold_xval';
else
    S.thisSelector = 'TrainTestOneIterGroup';
    S.thisSigIntenseSelector = 'TrainTestOneIterGroup';
end

%% Smoothing Parameters
S.use_unsmoothed = 1;
S.smoothTxt = {'smoothed' 'unsmoothed'};

if S.use_unsmoothed
    par.filesForPatterns = par.wascanfiles;
else
    par.filesForPatterns = par.swascanfiles;
end

S.filenames_train = vertcat(par.filesForPatterns{S.TrainRuns})  ;
S.filenames_test = vertcat(par.filesForPatterns{S.TestRuns}) ;

S.filenames = vertcat(S.filenames_train, S.filenames_test);

S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters
S.runs_vector =  [par.numvols(S.TrainRuns) par.numvols(S.TestRuns)];
S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = 2;

%% Volume Parameters
S.vol_info = spm_vol(fullfile(S.univar_dir, 'beta_0001.hdr')); %get functional data resolution info for spm .img writing
S.roi_name = 'tsrmask23.hdr';
S.roi_file = ['/Users/gordonam/Studies/AG1/mvpa/mask_012510/' S.roi_name];

S.secondaryMask = [];

%% Feature Selection Parameters
S.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)

%% Workspace Parameters
S.use_premade_workspace = 1;
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '_train_' S.trainTask '_test_' S.testTask '.mat']);
S.preprocPatName = 'spiral_d_z';
S.preprocPatCondensedName = 'spiral_d_z_condensed';

if isempty(S.secondaryMask)
    S.preprocPatNameMasked = 'spiral_d_z';
else
    S.preprocPatNameMasked = 'spiral_d_z_masked';
end

%% Mean Signal Extraction Params
S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];

S.ROI1_file = '/Users/gordonam/Studies/AG1/mvpa/VOTC_masks/rnondialtedHipp.m.hdr';
S.ROI2_file = '/Users/gordonam/Studies/AG1/mvpa/VOTC_masks/rnondialtedHipp.m.hdr';

S.logreg_2Features = 0;
S.defineROIsFromANOVAFS = 0;

S.zscoreIntensityVals = 0;
%% Iteration Parameters
S.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 1; % equate number of trials in conditions 1 and 2
S.numBalancedIts = 1;

%% Z-Scoring and outlier detection
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.importance_maps_dir=[S.expt_dir 'mvpa_results/ImpMaps_corr_hipp' date '/' S.roi_name;  ];
S.generate_importance_maps = 0;
S.generateBetaMaps = 0;

%% Special types of analysis
S.betaSeriesAnalysis = 0;
S.searchlightAnalysis = 0;
S.extractMeanSignal = 1;
S.nMRMR=0;
S.nHierClusts=0;

S.connectivityInteraction = 0;
S.ROIconnectivity =  '/Users/gordonam/Studies/AG1/mvpa/VOTC_masks/rnondialtedHipp.m.hdr';
S.corrClassification = 0;

%% TR Weighting

%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier
S.TR_weights_set = {{[0 0 0 .5 .5 ] [0 0 0 .5 .5 ]}};

%% classifier parameters
S.classify = 1;
S.class_args.train_funct_name = 'train_liblinear';
S.class_args.test_funct_name = 'test_liblinear';
S.class_args.nHidden = 0;
S.class_args.pThresh = 1;
S.class_args.constant = true;
S.class_args.prefitWeights = true;
S.perfmet_functs =  'perfmet_maxclass';
S.statmap_funct =  'AG_statmap_anova';

S.pca_proportion_var = 0; %0 for no pca on classification data.  otherwise, pca will be run.  specify what proportion of the data you want explained (e.g. for .99, the first n components will be used, where 99% of the variance is described by components 1 through n).
S.pca_nComps = 0;
S.nVoxSet = 0;  %S.Ps(pP); % 0 = no feature selection

S.class_args.libLin = '-q -s 0 -B 1';
S.class_args.libsvm = '-q -s 0 -t 1 -d 2';
S.class_args.chooseOptimalPenalty = 0;
S.class_args.penaltyRange = [.1];
S.class_args.nFolds = 10;
S.class_args.penalty = .1;
S.class_args.classType = 'libLin';
S.class_args.radialBasisSelection = 0;
S.class_args.nFoldsPenaltySelection = 10;

S.penaltyParams = .1;
S.nPlsCompsSet = 0;
S.class_args.nVox = 0;
S.linReg = 0;

%% pattern correlation
S.correlatePatterns = 0;
end