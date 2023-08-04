
%% model

% Assuming your data is in a table named 'data' with columns 'pupilSize', 'animalID', and 'trialNumber'

% Calculate the average pupil size for each trial and animal
averagePupilSize = grpstats(data, {'animalID', 'trialNumber'}, 'mean', 'DataVars', 'pupilSize');

% Group the data by trial and calculate the mean pupil size
meanPupilSizeByTrial = grpstats(averagePupilSize, 'trialNumber', 'mean', 'DataVars', 'mean_pupilSize');

% Create a matrix to store the average pupil size for each trial
averagePupilMatrix = reshape(meanPupilSizeByTrial.mean_pupilSize, [], 12);

% Perform linear regression for each direction separately
slopeCoefficients = zeros(1, 12);
pValues = zeros(1, 12);

for i = 1:12
    model = fitlm(meanPupilSizeByTrial, ['mean_pupilSize_' num2str(i) ' ~ trialNumber']);
    slopeCoefficients(i) = model.Coefficients.Estimate(2);
    pValues(i) = model.Coefficients.pValue(2);
end

% Display the results
disp('Slope coefficients:');
disp(slopeCoefficients);

disp('p-values:');
disp(pValues);