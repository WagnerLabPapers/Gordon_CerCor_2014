function gpar = AG1GroupParams(subjArray)
% parameters for group univariate analyses
% amg

gpar.subjArray = subjArray;
gpar.conGroup = { 'RetAnalysis' };
gpar.tasks = {'Ret'};
gpar.expt_dir = '/Users/gordonam/Studies/AG1/fmri_data2/';
gpar.modelTemplate = '/Users/gordonam/Studies/AG1/scripts/GroupTemplate8.mat';
gpar.constat = 'T';
gpar.exMask = [];

for t = 1:length(gpar.tasks)

    % use the first subject as a template
    gpar.task{t}.conTemplate =     fullfile(gpar.expt_dir, 'ag1_021509', gpar.conGroup{t}, 'SPM.mat');
    ldTemp = load(gpar.task{t}.conTemplate);
    gpar.task{t}.SPMcons = ldTemp.SPM.xCon;
    
    % for each contrast of interest, shuffle in contrasts from each subject
    % in subjArray
    for c= 1:length(gpar.task{t}.SPMcons);
        
        % contrast directory and name
        gpar.task{t}.cons{c}.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t},gpar.task{t}.SPMcons(c).name)};
        gpar.task{t}.cons{c}.name = gpar.task{t}.SPMcons(c).name;
        
        % put each subject's contrast image in the gpar structure.
        for s = 1:length(gpar.subjArray)
            gpar.task{t}.cons{c}.scans{s} = fullfile(gpar.expt_dir, gpar.subjArray{s}, gpar.conGroup{t}, ['con_' prepend(num2str(c), 4) '.img']);
        end
        
    end
end