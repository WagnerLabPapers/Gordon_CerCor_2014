function AG1setcontrasts(subpar, task)

% set contrasts for a single subject spm contrast
% written by jbh
% modified by amg



% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end

cd(par.(task).dir);
fprintf('\nLoading SPM...');
load SPM
fprintf('Done');


% T-contrasts
%---------------------------------------------------------------------------

if strcmp(task, 'Enc')
    cnames  = {'PConfVsSConf'; 'SConfVsPConf'; 'PAllVsSAll'; 'SAllVsPAll'; 'ConfVsNConf';...
        'NConfVsConf'; 'PConfVsPNConf'; 'PNConfVsPConf'; 'SConfVsSNConf'; 'SNConfVsSConf';...
        'PNConfVsSNConf'; 'SNConfVsPNConf'};
    
    cvals = {[1 0 0 -1  0 0 0 0 0 0 0 0]; [-1 0 0 1 0 0 0 0 0 0 0 0]; [1 0 0 -1 0 0 1 0 0 -1 0 0];...
        [-1 0 0 1 0 0 -1 0 0 1 0 0 ]; [1 0 0 1 0 0 -1 0 0 -1 0 0]; [-1 0 0 -1 0 0 1 0 0 1 0 0];...
        [1 0 0 0 0 0 -1 0 0 0 0 0]; [-1 0 0 0 0 0 1 0 0 0 0 0]; [0 0 0 1 0 0 0 0 0 -1 0 0]; ...
        [0 0 0 -1 0 0 0 0 0 1 0 0]; [0 0 0 0 0 0 1 0 0 -1 0 0]; [0 0 0 0 0 0 -1 0 0 1 0 0]};
elseif strcmp(task, 'Ret')
    cnames = {'Cor_Vs_Inc' 'P_Vs_S' 'Mem_X_Task' 'PCor_Vs_SCor'...
        'PInc_Vs_SInc' 'PCor_Vs_PInc' 'SCor_Vs_SInc'...
        'All_Vs_Fix' 'AllCor_Vs_Fix' 'P_Vs_Fix' 'S_Vs_Fix' 'PCor_Vs_Fix'...
        'SCor_Vs_Fix'};
    
    cvals = {[1 1 -1 -1 0 ] [1 -1 1 -1 0] [1 -1 -1 1 0] [1 -1 0 0] ...
        [0 0 1 -1] [1 0 -1 0 0] [0 1 0 -1 0]...
        [1 1 1 1 1] [1 1 0 0 0] [1 0 1 0 0] [0 1 0 1 0] [1 0 0 0 0] [0 1 0 0 0]};
end

% preallocate
con_name(1:length(cnames)) = {''};
con_vals = cell(1, length(cnames));

for Tt = 1:length(cnames)
    
    % make names
    con_name{Tt} = cnames{Tt};
    
    %puts two zeros between each element of cvals (to account for
    %inclusion of time and dispersion derivatives).
    val_processor_h = vertcat(cvals{Tt}, zeros(2, length(cvals{Tt})));
    val_processor = horzcat(val_processor_h(:))';
    
    
    con_vals{Tt} = val_processor;
end


% put contrasts into SPM/write to file
fprintf('\nBeginning contrasts on subject %s\n', par.substr);


cfid = fopen('conlist','wt');
fprintf(cfid, 'Contrasts for Sub %s\nLast run on %s\n', par.substr, date);

% Loop over created contrasts
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
copyfile('conlist',[par.logdir filesep 'conlist-' date]);

% Change back directory
cd(origdir);
return;

