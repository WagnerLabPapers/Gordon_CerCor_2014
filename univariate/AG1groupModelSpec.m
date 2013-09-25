function AG1groupModelSpec(par)
%
% Sets up the group-level spm design matrix and runs a design.
% 
% written by amg

for t = 1:length(par.tasks)
    clear jobs
    
    % load the template
    load (par.modelTemplate)
    yTemplate = load(fullfile(par.expt_dir, 'ag1_021509', par.conGroup{t}, 'SPM'));
    yCons = yTemplate.SPM.xCon;
    
    % loop through each of the contrasts
    for c= 1:length(yCons);
        
        jobs{c}.stats{1}.factorial_design.des.t1.scans = [];
        
        % make a directory for each contrast
        if ~exist(fullfile(par.expt_dir, 'group_analyses', par.tasks{t}, par.task{t}.SPMcons(c).name));
            mkdir(fullfile(par.expt_dir, 'group_analyses', par.tasks{t}, par.task{t}.SPMcons(c).name));
        end
        
        % shuffle in design stuff from group param struct
        jobs{c}.stats{1}.factorial_design.dir = {fullfile(par.expt_dir, 'group_analyses', par.tasks{t}, par.task{t}.SPMcons(c).name)};
        jobs{c}.stats{1}.factorial_design.masking.em{1} = par.exMask;
        jobs{c}.stats{1}.factorial_design.cov(1).cname = par.covName;
        jobs{c}.stats{1}.factorial_design.cov(1).c = par.covVec;
        jobs{c}.stats{1}.factorial_design.cov(1).iCC = 1;
        jobs{c}.stats{1}.factorial_design.cov(1).iCFI = 1;
        
        % place each subject's contrast image into the jobs structure.
        for s = 1:length(par.subjArray)
            jobs{c}.stats{1}.factorial_design.des.t1.scans{s} = fullfile(par.expt_dir, par.subjArray{s}, par.conGroup{t}, ['con_' prepend(num2str(c), 4) '.img']);
        end
    end
    
    % run the job for this task.
    spm_jobman('run',jobs)
    
end





        
        
        
        
        
