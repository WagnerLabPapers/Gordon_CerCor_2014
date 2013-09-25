function [groupDat groupIdx] = AG1GroupBehAnalysis(subjArray)

% summarize behavioral data across all subjects
% written by amg

for s = 1:length(subjArray)
    subject = subjArray{s};
    [dat idx]  = AG1SingleSubjBehAnalysis(subject);
    
    groupDat{s} = dat;
    groupIdx{s} = idx;
    
end