function [res, results]= AG_run_mvpa_general(subj_array, saveName)

% run an MVPA analysis 
% written by amg
% modified from code by jr


%%
for b=(1:length(subj_array))
    
    %% setup
    S = AG_mvpa_params(subj_array{b}); % load params
    tic; % start timing
    S.saveName= saveName; %name of classifier results
    S.subj_array = subj_array; %subjects
   
    %% Which TRs will be subjected to classification?
    S.TR_weights_train = S.TR_weights_set{1};
    S.TR_weights_test = S.TR_weights_set{2};
    
    S.TRs_to_average_over = 1:length(S.TR_weights_test);
    S.TRs_to_average_over_train = 1:length(S.TR_weights_train);
    S.TRs_to_average_over_test = 1:length(S.TR_weights_test); 
    
    
    %% Onsets
    S.thisTrain = S.trainTask;
    S.thisTest = S.testTask;
    S.univar_dir.enc = [S.expt_dir '/' S.subj_id '/encAnalysis'];
    S.univar_dir.ret = [S.expt_dir '/' S.subj_id '/retAnalysis'];
    S.idxThisTrain = 1;
    S.idxThisTest = 1;
    [S.onsets, S.img_files] = AG_mvpa_onsets_and_images(S);
    S.num_conds = size(S.onsets,2);
    
    %% Workspace stuff - load a premade workspace to save on time
    existWorkspace = exist(S.workspace);
    
    % load workspace
    if (S.use_premade_workspace&&existWorkspace)
        load(S.workspace, 'subj');
    else
        [subj] = PM_mvpa_load_and_preprocess_raw_data(S);
    end
    
    %% mask a workspace mask with another mask.
    if ~isempty(S.secondaryMask)
        subj = load_spm_mask(subj, 'secondaryMask', S.secondaryMask);
        subj = intersect_masks(subj,S.roi_name,'secondaryMask');
        subj = create_pattern_from_mask(subj, 'spiral_d_z', subj.masks{end}.name , 'spiral_d_z_masked');
    end
    
    %% for each iteration of a full classification
    for n = 1: S.num_results_iter
        
         % initialize regs matrix as conditions x timepoints
        all_regs = zeros(S.num_conds,S.num_vols);
        
        % convert from seconds to TRs
        for cond = 1: S.num_conds
            for trial = 1: length(S.onsets{cond})
                time_idx = round(S.onsets{cond}(trial)/S.TR) + 1; 
                all_regs(cond, round(time_idx)) = 1;
            end
        end
        
        % condense regs by removing zeros
        condensed_runs = [];
        condensed_regs_of_interest = [];
        trial_counter = 1;
        for i = 1: size(all_regs,2)
            if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
                condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
                condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                trial_counter = trial_counter + 1;
            end
        end
        idx_condense =find(sum(all_regs));
        
        % condense across runs
        trial_idx = 1;
        m_runs = 0;
        for r = 1:length(S.meta_runs)
            m_runs(trial_idx:trial_idx+S.meta_runs(r)-1)=r;
            trial_idx = trial_idx+S.meta_runs(r);
        end
        meta_runs_condensed = m_runs(idx_condense);
        
        %Train on one set of the data, test on the other
        TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.TrainRuns)) + ...
            2*ismember(meta_runs_condensed, length(S.TrainRuns)+1: length(S.runs_vector));
        
        % which trials to include
        actives_h = ones(size(meta_runs_condensed));
        actives = actives_h;
        subj = init_object(subj,'selector','actives');
        subj = set_mat(subj,'selector','actives', actives);
        
        all_trials = sum(all_regs,1);
        

        
        if strcmp(S.trainTask, S.testTask)
            % applies when train and test data are the same
            meta_runs_train = find(all_trials);
            meta_runs_test = [];
        else
            %applies only when train and test data are different
            meta_runs_train = idx_condense(find(TrainTestOneIter==1));
            meta_runs_test = idx_condense(find(TrainTestOneIter==2));
        end
        
        %% load data
        subj.patterns{1}.mat = [];

        % training set data
        data_by_TR_train = [];
        for dt = 1:length(S.TR_weights_train)
            data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:,meta_runs_train+(dt-1));
        end
        temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
        clear data_by_TR_train
        
        % testing set data
        data_by_TR_test = [];
        for dt = 1:length(S.TR_weights_test)
            data_by_TR_test(dt,:,:) = S.TR_weights_test(dt)*subj.patterns{end}.mat(:,meta_runs_test+(dt-1));
        end
        temporally_condensed_data_test = squeeze(sum(data_by_TR_test(S.TRs_to_average_over_test,:,:),1));
        clear data_by_TR_test
        
        % concatenate training and testing data
        temporally_condensed_data = horzcat(temporally_condensed_data_train, temporally_condensed_data_test);
        clear temporally_condensed_data_train;
        clear temporally_condensed_data_test;
        
        %% modify selectors in subj object
        
        %create TrainTestOneIter selector
        subj = init_object(subj,'selector','TrainTestOneIter');
        subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
        subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
        
        %create meta_runs_condensed selector
        subj = init_object(subj,'selector','meta_runs_condensed');
        subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
        subj = create_xvalid_indices(subj,'meta_runs_condensed');
        
        % create random selector, with S.nFolds iterations
        randomNFold = ceil(shuffle(1:length(meta_runs_condensed))/(length(meta_runs_condensed)/S.nFolds));
        subj = init_object(subj,'selector','randomNFold');
        subj = set_mat(subj,'selector','randomNFold', randomNFold);
        subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
        subj = PM_create_xvalid_indices_trainActivesOnly(subj,'randomNFold');
        
        % add S.thisSelector to subj object
        grp = find_group(subj, 'selector', S.thisSelector);
        for g = 1:length(grp)
            this_mat = get_mat(subj,'selector',grp{g});
            this_mat(this_mat==1) = this_mat(this_mat==1) .* actives(this_mat==1);
            subj = set_mat(subj,'selector',grp{g},this_mat);
        end
        
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);
        S.classSelector = S.thisSelector;
        
        %% regressors in subj object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
        subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
        
        %% patterns in subj object
        if ~isempty(S.secondaryMask)
            subj = duplicate_object(subj,'pattern','spiral_d_z_masked',S.preprocPatCondensedName);
        else
            subj = duplicate_object(subj,'pattern',S.preprocPatName,S.preprocPatCondensedName);
        end
        subj.patterns{end-1}.mat = [];
        subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
        zhist = sprintf('Pattern ''%s'' created by AG custom code',S.preprocPatCondensedName);
        subj = add_history(subj,'pattern',S.preprocPatCondensedName,zhist,true);
        
        % clean up workspace to save RAM
        subj = remove_mat(subj,'pattern',S.preprocPatName);
        
        %% equate the number of stimuli in each class of the training set.
        if S.equate_number_of_trials_in_cond_1_and_2
            subj = PM_balanceTrainPats(S, subj);
            S.classSelector = [S.thisSelector 'balanced'];
        end
                
        %% run feature selection ANOVA: specify #of voxels (if desired)
        if S.class_args.nVox>0
            subj = AG_feature_select_top_N_vox(subj,S.preprocPatCondensedName,'conds',S.classSelector, ...
                'nVox_thresh',S.class_args.nVox, 'statmap_funct', S.statmap_funct, 'statmap_arg',statmap_arg);
            
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        end
        classifier_pattern = S.preprocPatCondensedName;
        
        %% begin the classification
        S.class_args.penalty = S.penaltyParams;
        if S.classify
            if S.extractMeanSignal
                % instead of doing classification, just extract
                % the mean signal in a given set of regions and
                % compare them
                [subj results] =  AG_extractMeanSignal(subj,S);
            else
                %run the classification.
                [subj results] = cross_validation(subj,classifier_pattern,'conds', ...
                    S.classSelector, classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
            end
            
            %set up importance maps.
            if S.generate_importance_maps == 1
                for rif = 1:length(results.iterations);
                    results_IW{rif}.iterations(1).scratchpad.net.IW{1} = ...
                        results.iterations(rif).scratchpad.model.w(1:(end-1))';
                end
            end
            
            %store results
            res.subj{b}.iter{n} = results;
        end

        
        %% store and save things
                
        res.subj{b}.S = S;
        res.subjArray = subj_array;
        
        if ~(exist(S.group_mvpa_dir))
            mkdir(S.group_mvpa_dir);
        end
        
        % save the results structure
        save (fullfile(S.group_mvpa_dir, S.saveName), 'res');
        
        time2finish = toc/60;
        display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
        
        %% generate importance maps.
        if S.generate_importance_maps
            PM_generate_importance_maps(subj, results, results_IW, S)
        end
    end
end


















