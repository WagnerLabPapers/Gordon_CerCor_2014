function [subj] = AG1mvpa_load_and_preprocess_raw_data(S)

% written by amg

% create runs index
trial_idx = 1;
for r = 1:length(S.runs_vector)
    runs(trial_idx:trial_idx+S.runs_vector(r)-1)=r;
    trial_idx = trial_idx+S.runs_vector(r);
end
S.runs = runs;

% load patterns
subj = init_subj(S.exp_name,S.subj_id);
subj = load_spm_mask(subj,S.roi_name,S.roi_file);
subj = load_analyze_pattern(subj,'spiral',S.roi_name, S.img_files,'single',true);

% set runs
subj = init_object(subj,'selector','runs');
subj = set_mat(subj,'selector','runs',runs);

% detrend runs
subj = detrend_runs(subj,'spiral','runs'); 
subj = remove_mat(subj,'pattern','spiral');

% zscore the data
subj = zscore_runs(subj,'spiral_d','runs');
subj = remove_mat(subj,'pattern','spiral_d');

%save the workspace
if ~exist(S.workspace_dir)
    mkdir(S.workspace_dir);
end

cd (S.workspace_dir);
save (S.workspace);

end


