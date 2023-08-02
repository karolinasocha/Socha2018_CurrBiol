
%% select one session with weight probability
dataArray = [10, 20, 30, 40]; % Example array
weights=(1/length(dataArray))*(ones(1,numel(dataArray)))

numSamples = 1; % number of session selected
selectedIndices = randsample(numel(dataArray), numSamples, true, weights);

datasample