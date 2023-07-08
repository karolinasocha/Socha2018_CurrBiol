%% loading data from file
% loading data data_behav_40expt 40 sessions for behavioral analysis only
% protocol: not-randomized stimulus

% 08-07-2023: added animal_id and session_id:
% - to be save together as save_behav_40expt 
% - for future statistical tests
% 
% data necessary to generate plots for Figure1

clear all

% changes made 08-07-2023 by KSocha
% this is old directory: 
% load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')
% new directory: 
% load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;
diam=diameter;

%% added number of animals
animalList={}
for iAn=1:size(data_behav_40expt.newdir,1)
    animalList{iAn}=data_behav_40expt.newdir{iAn}(8:12);
end

animalUnique=unique(animalList)
j=0;
for iAn=1:size(animalUnique,2)
    j=j+1
    indices = find(strcmp(animalList, animalUnique{iAn}));
    for ii=1:length(indices);
    data_behav_40expt.animal_id(indices(ii))=j;
    end   
end

%% added sessoin id to each animal

for iAn=1:size(animalUnique,2)
    indices = find(strcmp(animalList, animalUnique{iAn}));
    j=0
        for ii=1:length(indices);
            j=j+1
            data_behav_40expt.session_id(indices(ii))=j;
        end   
end

save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\data_behav_40expt_v2.mat','data_behav_40expt')
 
