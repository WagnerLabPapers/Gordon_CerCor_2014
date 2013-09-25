function [ dat ] = ClassifierPermutationTest( res )
%Given a classifier output strucutre, perform a permutation test to
%determine the significance of the classification performance

% written by amg

nIts = 1000;   %number of iterations
baseLine = .5;  %assumed base rate of accuracy
condsOfInterest = {'PersonConf' 'SceneConf'}; % conditions of interest

for s = 1:length(res.subj)
    fprintf('running subject %g of %g \n', s, length(res.subj));
    
    Sres = res.subj{s}.penalty.nVox.weights.S;
    
    %% find the onsets you're interested in;
    mvpa_ons = load(fullfile(Sres.onsetsTestDir, 'mvpa_ons'));
    ixCondsOfInterest = ismember(mvpa_ons.names, condsOfInterest);
    
    %% store relevant onsets
    if Sres.sameTaskTrainAndTest
        relevantOnsets = vertcat(mvpa_ons.onsets{ixCondsOfInterest});
    else
        relevantOnsets = Sres.durTrain + vertcat(mvpa_ons.onsets{ixCondsOfInterest});
    end
    classifierOnsets = sort(vertcat(Sres.onsets{:}));
    testOnsets = classifierOnsets([res.subj{s}.penalty.nVox.weights.iter{1}.iterations(:).test_idx]);
    ixRelevantOnsets = ismember(testOnsets, relevantOnsets)';
    
    %% extract guesses and actual labels of the given classification
    resBase = res.subj{s}.penalty.nVox.weights.iter{1}.iterations;
    for i=1:length(resBase)
        desireds_h{i} = resBase(i).perfmet.desireds;
        guesses_h{i} = resBase(i).perfmet.guesses;
    end
    desireds = [desireds_h{:}];
    guesses = [guesses_h{:}];
    
    % OI = values of interest
    desiredsOI = desireds(ixRelevantOnsets);
    guessesOI = guesses(ixRelevantOnsets);
    classPerfOI =  mean(desiredsOI==guessesOI);
    
    %% mean accuracy across all groups
    for j = unique(desiredsOI)
        classPerfOI_meanAcrossGroups_h(j) = mean(guessesOI(desiredsOI==j)==j);
    end
    classPerfOI_meanAcrossGroups = mean(classPerfOI_meanAcrossGroups_h);
    
    %% permute the guesses nIts times
    dat.sub(s).perfDist = nan(1,nIts);
    for n=1:nIts
        thisPerm = shuffle(guessesOI);
        dat.sub(s).perfDist(n) = mean(thisPerm==desiredsOI);
    end

    %% create a distribution of pvalues
    distPval = mean(classPerfOI<=dat.sub(s).perfDist);
    if distPval==0
        dat.sub(s).pVal = 1/nIts;
    else
        dat.sub(s).pVal = distPval;
    end
    
    %% bootstrapped statistic
    for j=unique(desiredsOI)
        classMembers = (desiredsOI==j);               
        bootstat_h = bootstrp(nIts, (@(x) x), guessesOI);
        bootstat_h2(:,j) = mean(bootstat_h(:,classMembers)==j,2);
    end
    bootstat = mean(bootstat_h2,2);
    
    %% put relevant variables in the dat structure
    dat.sub(s).perfByClass = classPerfOI_meanAcrossGroups_h;
    dat.sub(s).actualPerf = classPerfOI_meanAcrossGroups;
    dat.sub(s).empiricalMedianPerm = median(dat.sub(s).perfDist);
    dat.sub(s).empiricalMedianBoot = median(bootstat);
    dat.sub(s).diffStatPerm = classPerfOI_meanAcrossGroups - median(dat.sub(s).perfDist);
    dat.sub(s).diffStatBoot= classPerfOI_meanAcrossGroups - median(bootstat);
end

%% group level statistics
[~, dat.group.tPermute.p, ~, dat.group.tPermute.stat]  = ttest([dat.sub(:).diffStatPerm]);
[~, dat.group.tassumeBaseline.p, ~, dat.group.tassumeBaseline.stat]  = ttest([dat.sub(:).actualPerf], baseLine);
[~, dat.group.tBoot.p, ~, dat.group.tBoot.stat]  = ttest([dat.sub(:).diffStatBoot]);

