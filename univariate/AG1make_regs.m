function dat = AG1make_regs(par, saveit)

% makes regressors relevant to experiment
% written by amg


%% general information 
if ~exist('saveit')
    saveit = input('\n save onsets and regressors?  1 = yes  0 = no \n \n');
end

expDir = par.fmridir;

myDir.subj = fullfile(expDir, par.substr);
myDir.Beh= (fullfile(expDir, par.substr, 'behav'));

cd (myDir.Beh);

%data file names
encDirH = dir ('AG1_enc*');
retDirH = dir ('AG1_ret*');

%load the raw data files
rawDat.Enc = load(encDirH.name);
rawDat.Ret = load(retDirH.name);

[dat idx] = AG1SingleSubjBehAnalysis;


%% Encoding regressors
clear ons onsets durations names;

ESN = length(par.Enc.Scans);

onsh.Enc.alltrials = [];
onsh.Enc.h = rawDat.Enc.EncData(ESN).onset;

% trial onsets
for i = 1:ESN
    onsh.Enc.alltrials = round([onsh.Enc.alltrials, (onsh.Enc.h-11 +par.Enc.dur*(i-1) )]);
end

% session regressors
for i = 1:ESN-1
    sessReg.Enc((i-1)*par.Enc.nTRs+1:i*par.Enc.nTRs,i) = ones(par.Enc.nTRs,1);
end

% onsets of interest
ons.Enc.personConf = onsh.Enc.alltrials(find(idx.enc.personConf));
ons.Enc.sceneConf = onsh.Enc.alltrials(find(idx.enc.sceneConf));
ons.Enc.junk = onsh.Enc.alltrials(find(idx.enc.junk));

% remove junk onset if there are no junk trials
if(isempty(ons.Enc.junk))
    ons.Enc = rmfield(ons.Enc, 'junk');
end

% reshape onsets, durations, and names
fn = fieldnames(ons.Enc);
allEncOns = [];
for f = 1:length(fn);
    if isempty(getfield(ons.Enc, fn{f}))
        fprintf('%s has empty Encoding condition %s', par.substr, fn{f});
    end
    onsets{f} = ons.Enc.(fn{f});
    durations{f} = 0;
    names{f} = fn{f};
    allEncOns = [allEncOns, onsets{f}];
end

if length(allEncOns)~=par.Enc.nTrials
    fprintf( '\n Warning: %s encoding onsets for subject %s \n', num2str(length(allEncOns)), par.substr);
end

%% Encoding covariates

% sessions
sessReg.Enc = zeros(par.Enc.nTRs*ESN,ESN-1);
for i = 1:ESN-1
    sessReg.Enc((i-1)*par.Enc.nTRs+1:i*par.Enc.nTRs,i) = ones(par.Enc.nTRs,1);
end

% artifacts
Arts.Raw = load(fullfile(par.artrepdir, ['art_global_modified_' par.substr ]));
Arts.Enc.Raw = Arts.Raw.allArt(par.EncScans);
Arts.Enc.Mat = vertcat(Arts.Enc.Raw{:});

Arts.Enc.Idx = find(Arts.Enc.Mat==0);
Arts.Enc.Reg = zeros(length(Arts.Enc.Mat), length(Arts.Enc.Idx));
for i = 1:length(Arts.Enc.Idx)
    Arts.Enc.Reg(Arts.Enc.Idx(i),i) = 1;
end

% combine sessions and artifacts covariates
R = horzcat(sessReg.Enc, Arts.Enc.Reg);


ons.Enc.dur = 0;
dat.sub.Enc.ons = ons.Enc;
dat.sub.Enc.sessReg = sessReg.Enc;

% save the regressors and onsets in SPM format
myDir.Enc = (par.Enc.dir);
if saveit
    if ~exist(myDir.Enc)
        mkdir(myDir.Enc);
    end
    cd (myDir.Enc);
    
    save ons onsets durations names;
    save regs.mat R
    
end


%% Retrieval regressors
clear ons onsets durations names

% onsets for all trials
onsh.Ret.alltrials = [];
onsh.Ret.h = rawDat.Ret.retData(par.Ret.nSess).onset;
for i = 1:par.Ret.nSess
    onsh.Ret.alltrials = round([onsh.Ret.alltrials, (onsh.Ret.h-12 +par.Ret.dur*(i-1) )]);
end

% onsets for regressors of interest
ons.Ret.PCorConf = onsh.Ret.alltrials(find(idx.ret.PersonCorConf));
ons.Ret.SCorConf = onsh.Ret.alltrials(find(idx.ret.SceneCorConf));
ons.Ret.PIncConf = onsh.Ret.alltrials(find(idx.ret.PersonIncConf));
ons.Ret.SIncConf = onsh.Ret.alltrials(find(idx.ret.SceneIncConf));

% all trials that aren't in a regressor of interest are junk trials
idx.ret.Junk = ones(par.Ret.nTrials,1) - idx.ret.PersonCorConf - idx.ret.SceneCorConf...
    -idx.ret.PersonIncConf - idx.ret.SceneIncConf;
ons.Ret.Junk = onsh.Ret.alltrials(find(idx.ret.Junk));

% shuffle onsets, durations, and names into a format SPM can read
allRetOns = [];
fn = fieldnames(ons.Ret);
for f = 1:length(fn);
    if isempty(getfield(ons.Ret, fn{f}))
        fprintf('%s has empty Retrieval condition %s', par.substr, fn{f});
    end
    onsets{f} = ons.Ret.(fn{f});
    durations{f} = 0;
    names{f} = fn{f};
    allRetOns = [allRetOns, onsets{f}];
end

% if the number of trials is different from the standard, throw a warning
if length(allRetOns)~=par.Enc.nTrials
    fprintf( '\n Warning: %s retrieval onsets for subject %s \n', num2str(length(allEncOns)), par.substr);
end

%% Retrieval covariates

% session covariate
sessReg.Ret = zeros((par.Ret.nSess-1)*par.Ret.nTRs,(par.Ret.nSess-1));
for i = 1:(par.Ret.nSess-1)
    sessReg.Ret((i-1)*par.Ret.nTRs+1:i*par.Ret.nTRs,i) = ones(par.Ret.nTRs,1);
end
dat.sub.Ret.ons = ons.Ret;
dat.sub.Ret.sessReg = sessReg.Ret;

% artifact covariate
Arts.Raw = load(fullfile(par.artrepdir, ['art_global_modified_' par.substr ]));
Arts.Ret.Raw = Arts.Raw.allArt(par.RetScans);
Arts.Ret.Mat = vertcat(Arts.Ret.Raw{:});
Arts.Ret.Idx = find(Arts.Ret.Mat==0);
Arts.Ret.Reg = zeros(length(Arts.Ret.Mat), length(Arts.Ret.Idx));
for i = 1:length(Arts.Ret.Idx)
    Arts.Ret.Reg(Arts.Ret.Idx(i),i) = 1;
end

% combine session and artifact covariates
R = horzcat(sessReg.Ret, Arts.Ret.Reg);
myDir.Ret = (fullfile(expDir, par.substr, 'retAnalysis'));

% save the onsets and regressors
if saveit
    if ~exist(myDir.Ret)
        mkdir(myDir.Ret);
    end
    cd (myDir.Ret);
    save ons.mat onsets durations names;
    save regs.mat R
end

    

