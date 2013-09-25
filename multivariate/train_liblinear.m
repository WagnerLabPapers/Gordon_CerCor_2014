
function [scratchpad] = train_liblinear(trainpats,traintargs,in_args)
%
%
% use liblinear toolbox to train a classifier with fmri data
% written by amg

v = (double(trainpats'));

for i=1:size(traintargs,2)
    choice(i,1) = find(traintargs(:,i))';
end


%%  pick the optimal cost parameter

choice_set = unique(choice);

if in_args.chooseOptimalPenalty
    
    
    idx.c1 = find(choice==choice_set(1));
    idx.c2 = find(choice==choice_set(2));
    
    guesses = nan(size(choice));
    
    %loop through indices, one at a time
    randomNFold1 = ceil(shuffle(1:length(idx.c1))/(length(idx.c1)/in_args.nFoldsPenaltySelection));
    randomNFold2 = ceil(shuffle(1:length(idx.c2))/(length(idx.c2)/in_args.nFoldsPenaltySelection));
    for s = 1:in_args.nFoldsPenaltySelection
        
       
        % omit trials to be tested on from training set.
        omitc1 = idx.c1(randomNFold1==s);
        omitc2 = idx.c2(randomNFold2==s);
        theseOmitted = [omitc1; omitc2];
       
        thisChoice = choice;
        thisChoice(theseOmitted) = [];
        
        thisV = v;
        thisV(theseOmitted,:) = [];
        
        
        % loop through all penalty parameters
        for i = 1:length(in_args.penaltyRange)
            l = in_args.penaltyRange(i);
            
            % liblinear training options.
            trainOpts_orig = in_args.libLin ;
            trainOptsOptimize = [trainOpts_orig ' -c ' num2str(l)];
            
            % train the classifier
            m = train(thisChoice, sparse(thisV), trainOptsOptimize);
            
            % test using held out (omitted) data
            [theseLabels,~]=predict(choice(theseOmitted), sparse(v(theseOmitted,:)), m);
            
        end
        % shuffle classifier guesses into this variable
        guesses(theseOmitted,i,r) = theseLabels;
    end
    
    
    choiceMat = repmat(choice,[1,length(in_args.penaltyRange), length(rSet)]);
    
    %% performance measures
    perf.tp = sum(guesses == choice_set(1) & choiceMat == choice_set(1));   % true pos
    perf.fp = sum(guesses == choice_set(1) & choiceMat == choice_set(2));  % false pos
    perf.fn = sum(guesses == choice_set(2) & choiceMat == choice_set(1));   % false neg
    perf.tn = sum(guesses == choice_set(2) & choiceMat == choice_set(2));  % true neg
    
    perf.Precision = perf.tp./(perf.tp+perf.fp);
    perf.Recall = perf.tp./(perf.tp+perf.fn);
    perf.TrueNegRate = perf.tn./(perf.tn+perf.fp);    
    perf.Accuracy = ((perf.tp)./(perf.tp + perf.fn) + (perf.tn)./(perf.tn + perf.fp))*.5;
    perf.F_Score = 2.*perf.Precision.*perf.Recall./ ...
        (perf.Precision+perf.Recall);
    
    % select the optimal parameter
    perfParam = perf.Accuracy;
    optPerf = max(max(perfParam));
    perfParam = squeeze(perfParam);
    idx.optPerf = find(perfParam==optPerf);
    opt_penalty = in_args.penaltyRange(idx.optPerf(1));

else
    % otherwise, don't optimize and just use what penalty parameter we've
    % specified
    opt_penalty = in_args.penalty;
end

scratchpad.opt_penalty = opt_penalty;
%% classify with the optimal penalty param, established in a non-biased
%% fashion by cross-validating the training data.

% train it, with leave one out cross validation
trainOpts_orig = in_args.libLin ;
trainOpts = [trainOpts_orig ' -c ' num2str(opt_penalty)];
model = train(choice, sparse(v), trainOpts);

% model information
scratchpad.classOrientation = model.Label';
scratchpad.logreg.betas = [model.w(end) model.w(1:(end-1))];
scratchpad.w(:,scratchpad.classOrientation(1),:) = scratchpad.logreg.betas';
scratchpad.w(:,scratchpad.classOrientation(2),:) = -1*scratchpad.logreg.betas';
scratchpad.classType = in_args.classType;
scratchpad.constant = in_args.constant;
scratchpad.choice = choice;
scratchpad.model = model;


