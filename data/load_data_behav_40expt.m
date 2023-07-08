%% loading data from file
% loading data data_behav_40expt 40 sessions for behavioral analysis only
% protocol: not-randomized stimulus

% data necessary to generate plots for Figure1

clear all

%load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

% 08-07-2023: added new directory to load data_behav_40expt
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;
diam=diameter;
session_id= data_behav_40expt.session_id;
animal_id= data_behav_40expt.animal_id;
