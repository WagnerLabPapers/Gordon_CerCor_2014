function [acts scratchpad] = test_liblinear(testpats,testtargs,scratchpad)

% use liblinear package to test classifier trained with fmri data with
% other fmri data
%
% written by amg

% should we add a constant term?
if scratchpad.constant
  testpats = [testpats; ones(1,cols(testpats))];
end

% sparsify the test patterns
thisTest = sparse(double(testpats));

% reshape the test targets data to fit liblinear format
for i=1:size(testTargsReoriented,2)
    theseTestLabels(i,1) = find(testtargs(:,i))';
end

% predict category, given test data.
[~, ~, probVals] = predict(theseTestLabels, thisTest', scratchpad.model, '-b 1');

% generate probabalistic classifier outcomes
acts = nan(size(probVals));
for i = 1:length(scratchpad.classOrientation)
    acts(:,i) = probVals(:,scratchpad.classOrientation(i),:);
end
acts = acts';
