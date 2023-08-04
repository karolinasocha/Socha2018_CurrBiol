% load pupil data zscored
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')
% pupil_data=[animal_id, sessions_id, trial_id, dir1,...,dir12, control_blank]

%% test across trials single animals
indecies=find(pupil_data.zscored_data(1,:)==1)
pupil_data_single=pupil_data.zscored_data(indecies,:)

numIterations=1000
nStim=12

bootstrap_results_trials = zeros(numIterations, 12);

for i = 1:numIterations
    bootstrap_sample = pupil_data(randi(numTrials, [numTrials, 1]), :);
    bootstrap_results_trials(i, :) = mean(bootstrap_sample);
end

% Calculate statistical significance and confidence interval for within-trial analysis
mean_trials = mean(pupil_data);
bootstrap_mean_trials = mean(bootstrap_results_trials);
bootstrap_std_trials = std(bootstrap_results_trials);
confidence_interval_trials = prctile(bootstrap_results_trials, [2.5, 97.5]);
p_value_trials = sum(bootstrap_results_trials >= repmat(mean_trials, [numIterations, 1])) / numIterations;

% Display the results
disp('Within-trial analysis:');
disp(['Mean pupil size: ', num2str(mean_trials)]);
disp(['Bootstrap mean pupil size: ', num2str(bootstrap_mean_trials)]);
disp(['Bootstrap standard deviation: ', num2str(bootstrap_std_trials)]);
disp(['Confidence interval (95%): ', num2str(confidence_interval_trials)]);
disp(['p-value: ', num2str(p_value_trials)]);

%% test across sessions

%% test across animals

% input: 
% pupil_size_data

%%

% Define the number of trials, sessions, and animals
numTrials = 10;
numSessions = % Specify the number of sessions in your data
numAnimals = % Specify the number of animals in your data

% Create the pupil_data matrix for within-trial analysis
pupil_data = % Initialize the matrix with the appropriate size

% Fill in the pupil_data matrix with your actual pupil size measurements
% Each row represents a trial, and each column represents a direction of gratings
% For example:
pupil_data = [
    % Trial 1
    pupil_size_1_direction_1, pupil_size_1_direction_2, ..., pupil_size_1_direction_12;
    % Trial 2
    pupil_size_2_direction_1, pupil_size_2_direction_2, ..., pupil_size_2_direction_12;
    ...
    % Trial numTrials
    pupil_size_numTrials_direction_1, pupil_size_numTrials_direction_2, ..., pupil_size_numTrials_direction_12
];

%%
% Define the number of bootstrap iterations
numIterations = 1000;

% Define the number of trials, sessions, and animals
numTrials = size(pupil_data, 1);
numSessions = size(pupil_data, 2) / 10; % Assuming 10 trials per session
numAnimals = size(pupil_data, 1) / numSessions; % Assuming equal sessions per animal

% Initialize matrices to store bootstrap results
bootstrap_results_trials = zeros(numIterations, 12);
bootstrap_results_sessions = zeros(numIterations, numSessions);
bootstrap_results_animals = zeros(numIterations, numAnimals);

% Perform bootstrapping for within-trial analysis
for i = 1:numIterations
    bootstrap_sample = pupil_data(randi(numTrials, [numTrials, 1]), :);
    bootstrap_results_trials(i, :) = mean(bootstrap_sample);
end

% Perform bootstrapping for within-session analysis
for i = 1:numIterations
    bootstrap_sample = pupil_data(randi(numSessions, [numSessions, 1]), :);
    bootstrap_results_sessions(i, :) = mean(bootstrap_sample);
end

% Perform bootstrapping for between-animal analysis
for i = 1:numIterations
    bootstrap_sample = pupil_data(randi(numAnimals, [numAnimals, 1]) + (0:numAnimals-1)*numSessions, :);
    bootstrap_results_animals(i, :) = mean(bootstrap_sample);
end

% Calculate statistical significance and confidence interval for within-trial analysis
mean_trials = mean(pupil_data);
bootstrap_mean_trials = mean(bootstrap_results_trials);
bootstrap_std_trials = std(bootstrap_results_trials);
confidence_interval_trials = prctile(bootstrap_results_trials, [2.5, 97.5]);
p_value_trials = sum(bootstrap_results_trials >= repmat(mean_trials, [numIterations, 1])) / numIterations;

% Calculate statistical significance and confidence interval for within-session analysis
mean_sessions = mean(reshape(pupil_data, numSessions, []));
bootstrap_mean_sessions = mean(bootstrap_results_sessions);
bootstrap_std_sessions = std(bootstrap_results_sessions);
confidence_interval_sessions = prctile(bootstrap_results_sessions, [2.5, 97.5]);
p_value_sessions = sum(bootstrap_results_sessions >= repmat(mean_sessions, [numIterations, 1])) / numIterations;

% Calculate statistical significance and confidence interval for between-animal analysis
mean_animals = mean(reshape(pupil_data, numSessions, numAnimals)');
bootstrap_mean_animals = mean(bootstrap_results_animals);
bootstrap_std_animals = std(bootstrap_results_animals);
confidence_interval_animals = prctile(bootstrap_results_animals, [2.5, 97.5]);
p_value_animals = sum(bootstrap_results_animals >= repmat(mean_animals, [numIterations, 1])) / numIterations;

% Display the results
disp('Within-trial analysis:');
disp(['Mean pupil size: ', num2str(mean_trials)]);
disp(['Bootstrap mean pupil size: ', num2str(bootstrap_mean_trials)]);
disp(['Bootstrap standard deviation: ', num2str(bootstrap_std_trials)]);
disp(['Confidence interval (95%): ', num2str(confidence_interval_trials)]);
disp(['p-value: ', num2str(p_value_trials)]);

disp('Within-session analysis:');
disp(['Mean pupil size: ', num2str(mean_sessions)]);
disp(['Bootstrap mean pupil size: ', num2str(bootstrap_mean_sessions)]);
disp(['Bootstrap standard deviation: ', num2str(bootstrap_std_sessions)]);
disp(['Confidence interval (95%): ', num2str(confidence_interval_sessions)]);
disp(['p-value: ', num2str(p_value_sessions)]);

disp('Between-animal analysis:');
disp(['Mean pupil size: ', num2str(mean_animals)]);
disp(['Bootstrap mean pupil size: ', num2str(bootstrap_mean_animals)]);
disp(['Bootstrap standard deviation: ', num2str(bootstrap_std_animals)]);
disp(['Confidence interval (95%): ', num2str(confidence_interval_animals)]);
disp(['p-value: ', num2str(p_value_animals)]);