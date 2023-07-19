
clear all
% zscored values takes average zscored per session
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')

session_list=unique(pupil_data.sessions_id)
animals_id=unique(pupil_data.animal_id)


%%
nconditions=13

data_all={}; % cells containing the parameters of inidividual areas

%data_in=pupil_average_data % averaged over sessions
%add_val=2;
data_in=pupil_data.zscored_data % all trials
add_val=3;    

for icond=1:nconditions
    data_all{icond}=data_in(:,icond+add_val);
end
%%
% Assuming you have a column vector 'data' with the measurements for each group
data = [data_all{1}; data_all{2}; data_all{3};...
    data_all{4} ; data_all{5}; data_all{6};...
    data_all{7};data_all{8};data_all{9};...
    data_all{10};data_all{11};data_all{12}]; % Replace with your actual data

% Create a grouping variable corresponding to each group
group = [repmat({'Group 1'}, numel(data_all{1}), 1);
         repmat({'Group 2'}, numel(data_all{2}), 1);
         repmat({'Group 3'}, numel(data_all{3}), 1);
         repmat({'Group 4'}, numel(data_all{4}), 1);
         repmat({'Group 5'}, numel(data_all{5}), 1);
         repmat({'Group 6'}, numel(data_all{6}), 1);
         repmat({'Group 7'}, numel(data_all{7}), 1);
         repmat({'Group 8'}, numel(data_all{8}), 1);
         repmat({'Group 9'}, numel(data_all{9}), 1);
         repmat({'Group 10'}, numel(data_all{10}), 1);
         repmat({'Group 11'}, numel(data_all{11}), 1);
         repmat({'Group 12'}, numel(data_all{12}), 1)]; % Replace with appropriate group labels

% Perform one-way ANOVA
[p, tbl, stats] = anova1(data, group);

% Display the ANOVA table and perform post-hoc tests if needed
disp(tbl);