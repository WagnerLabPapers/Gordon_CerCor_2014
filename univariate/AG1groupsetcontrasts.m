function AG1groupsetcontrasts(subpar, tsk, cnd)
% sets second-level spm contrasts
% written by amg
% based on code by jbh

origdir = pwd;


% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = AG1GroupParams(subpar);
end

STAT = par.constat;

cd(par.task{tsk}.cons{cnd}.dir{1})


fprintf('\nLoading SPM...');
load SPM
fprintf('Done');

% contrasts of interest
cnames = {par.task{tsk}.cons{cnd}.name; ['inv_' par.task{tsk}.cons{cnd}.name]};
cvals = {[1 0]; [-1 0]};

% preallocate
con_name(1:length(cnames)) = {''};
con_vals = cell(1, length(cnames));


for Tt = 1:length(cnames)
    % make names
    con_name{Tt} = cnames{Tt};
    con_vals{Tt} = cvals{Tt};
end


% put contrasts into SPM/write to file
fprintf('\nBeginning contrasts on task %s contrast %s\n', par.tasks{tsk}, par.task{tsk}.cons{cnd}.name);


cfid = fopen('conlist','wt');

% Loop over created contrasts
%-------------------------------------------------------------------
for k=1:length(con_vals)

    % Basic checking of contrast
    %-------------------------------------------------------------------
    [c,I,emsg,imsg] = spm_conman('ParseCon',con_vals{k},SPM.xX.xKXs,STAT);
    if ~isempty(emsg)
        disp(emsg);
        error('Error in contrast specification');
    else
        disp(imsg);
    end;

    % Fill-in the contrast structure
    %-------------------------------------------------------------------
    if all(I)
        DxCon = spm_FcUtil('Set',con_name{k},STAT,'c',c,SPM.xX.xKXs);
    else
        DxCon = [];
    end

    % Append to SPM.xCon. SPM will automatically save any contrasts that
    % evaluate successfully.
    %-------------------------------------------------------------------
    if isempty(SPM.xCon)
        SPM.xCon = DxCon;
    elseif ~isempty(DxCon)
        SPM.xCon(end+1) = DxCon;
    end
    SPM = spm_contrasts(SPM,length(SPM.xCon));
        
    fprintf(fopen('conlist','at'),'%d: %s\n%s\n\n',k, con_name{k},num2str(con_vals{k}));
end

fclose(cfid);


% Change back directory
cd(origdir);
return;


