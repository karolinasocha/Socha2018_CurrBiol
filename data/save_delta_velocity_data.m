%% create data set with nomalized to median pupil size
%% this is to set analysis on pupil size and assess an effect during stimulation

%% test between animals
%% add directory to functions
addpath(genpath('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\functions'))

%% Main figure 1
clear all
data_behav_40expt_directory='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\data_behav_40expt.mat'
load(data_behav_40expt_directory)
%load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;
diam=velocity;
session_id= data_behav_40expt.session_id;
animal_id= data_behav_40expt.animal_id;

%% creates epochs and re-order stimulus so it is consistent for all experiments 

for iAn=1:size(diameter,2)

diam_norm_median=(diam{iAn}-nanmedian(diam{iAn}))./nanmedian(diam{iAn});
diam_substracted_prctile=(diam{iAn}-prctile(diam{iAn},5));
diam_raw=diam{iAn};
diam_zscore=zscore(diam{iAn});

all_epochs=expt2{iAn}.frames.epochs;
stim_epochs=expt2{iAn}.frames.stims;

% raw data
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_raw;
[av_raw_diam_epochs{iAn},trials_raw_diam_epochs{iAn}, av_raw_diam_stims{iAn}, trials_raw_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_raw_diam_epochs_reordered{iAn}, trials_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_raw_diam_epochs{iAn},trials_raw_diam_stims{iAn});
[av_raw_diam_epochs_reordered{iAn}, av_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_raw_diam_epochs{iAn},av_raw_diam_stims{iAn});

% substracted from minimum value
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_substracted_prctile;
[av_prct_diam_epochs{iAn},trials_prct_diam_epochs{iAn}, av_prct_diam_stims{iAn}, trials_prct_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_prct_diam_epochs_reordered{iAn}, trials_prct_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_prct_diam_epochs{iAn},trials_prct_diam_stims{iAn});
[av_prct_diam_epochs_reordered{iAn}, av_prct_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_prct_diam_epochs{iAn},av_prct_diam_stims{iAn});


% norm to median
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_norm_median;
[av_median_diam_epochs{iAn},trials_median_diam_epochs{iAn}, av_median_diam_stims{iAn}, trials_median_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_median_diam_epochs_reordered{iAn}, trials_median_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_median_diam_epochs{iAn},trials_median_diam_stims{iAn});
[av_median_diam_epochs_reordered{iAn}, av_median_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_median_diam_epochs{iAn},av_median_diam_stims{iAn});

end

%% calculate average pupil size during epochs

for iAn=1:length(av_median_diam_epochs_reordered)
    
    [av_median_diam_stims_reordered_onset{iAn}, trials_median_diam_stims_reordered_onset{iAn},...
    av_median_diam_stims_reordered_offset{iAn}, trials_median_diam_stims_reordered_offset{iAn},...
    av_median_diam_stims_reordered_middle{iAn}, trials_median_diam_stims_reordered_middle{iAn},...
    av_median_diam_stims_reordered_entire{iAn}, trials_median_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_median_diam_stims_reordered{iAn}, trials_median_diam_stims_reordered{iAn} );

    [av_prct_diam_stims_reordered_onset{iAn}, trials_prct_diam_stims_reordered_onset{iAn},...
    av_prct_diam_stims_reordered_offset{iAn}, trials_prct_diam_stims_reordered_offset{iAn},...
    av_prct_diam_stims_reordered_middle{iAn}, trials_prct_diam_stims_reordered_middle{iAn},...
    av_prct_diam_stims_reordered_entire{iAn}, trials_prct_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_prct_diam_stims_reordered{iAn}, trials_prct_diam_stims_reordered{iAn} );

    [av_raw_diam_stims_reordered_onset{iAn}, trials_raw_diam_stims_reordered_onset{iAn},...
    av_raw_diam_stims_reordered_offset{iAn}, trials_raw_diam_stims_reordered_offset{iAn},...
    av_raw_diam_stims_reordered_middle{iAn}, trials_raw_diam_stims_reordered_middle{iAn},...
    av_raw_diam_stims_reordered_entire{iAn}, trials_raw_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_raw_diam_stims_reordered{iAn}, trials_raw_diam_stims_reordered{iAn} );

end

%%
for iAn=1:length(av_raw_diam_stims_reordered_entire)
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_raw_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_raw_diam_stims_reordered_offset{iAn};
    trials_raw_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    trials_raw_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    av_raw_diam_diff_relative_stimulation{iAn}=nanmean(trials_raw_diam_diff_relative_stimulation{iAn},1) ;
    
    
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_prct_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_prct_diam_stims_reordered_offset{iAn};
    trials_prct_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    %trials_prct_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_median_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_median_diam_stims_reordered_offset{iAn};
    trials_median_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    %trials_median_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    
end
% %% FIGURE
% for iAn=1:40
% figure(1)
% clf
% subplot(3,1,1)
% plot(squeeze(nanmean(trials_median_diam_difference_stimulation{iAn},1)))
% title('median norm delta')
% 
% subplot(3,1,2)
% plot(squeeze(nanmean(trials_prct_diam_difference_stimulation{iAn},1)))
% title('prct norm delta')
% 
% subplot(3,1,3)
% plot(squeeze(nanmean(trials_raw_diam_difference_stimulation{iAn},1)))
% title('raw delta')
% pause
% end
%%
clear new_velocity_data
new_velocity_data.newdir=data_behav_40expt.newdir;
new_velocity_data.stimulus=data_behav_40expt.stimulus;
new_velocity_data.velocity=data_behav_40expt.velocity;
new_velocity_data.expt=data_behav_40expt.expt;
new_velocity_data.expt2=data_behav_40expt.expt2;
new_velocity_data.diameter=data_behav_40expt.diameter;
new_velocity_data.diam=diameter;
new_velocity_data.sessions_id= data_behav_40expt.session_id;
new_velocity_data.animal_id= data_behav_40expt.animal_id;

new_velocity_data.av_median_diam_stims_reordered_onset=av_median_diam_stims_reordered_onset;
new_velocity_data.trials_median_diam_stims_reordered_onset=trials_median_diam_stims_reordered_onset;
new_velocity_data.av_median_diam_stims_reordered_offset=av_median_diam_stims_reordered_offset;
new_velocity_data.trials_median_diam_stims_reordered_offset=trials_median_diam_stims_reordered_offset;
new_velocity_data.av_median_diam_stims_reordered_middle=av_median_diam_stims_reordered_middle;
new_velocity_data.trials_median_diam_stims_reordered_middle=trials_median_diam_stims_reordered_middle;
new_velocity_data.av_median_diam_stims_reordered_entire=av_median_diam_stims_reordered_entire;
new_velocity_data.trials_median_diam_stims_reordered_entire=trials_median_diam_stims_reordered_entire;

new_velocity_data.trials_median_diam_difference_stimulation=trials_median_diam_difference_stimulation;

new_velocity_data.av_raw_diam_stims_reordered_onset=av_raw_diam_stims_reordered_onset;
new_velocity_data.trials_raw_diam_stims_reordered_onset=trials_raw_diam_stims_reordered_onset;
new_velocity_data.av_raw_diam_stims_reordered_offset=av_raw_diam_stims_reordered_offset;
new_velocity_data.trials_raw_diam_stims_reordered_offset=trials_raw_diam_stims_reordered_offset;
new_velocity_data.av_raw_diam_stims_reordered_middle=av_raw_diam_stims_reordered_middle;
new_velocity_data.trials_raw_diam_stims_reordered_middle=trials_raw_diam_stims_reordered_middle;
new_velocity_data.av_raw_diam_stims_reordered_entire=av_raw_diam_stims_reordered_entire;
new_velocity_data.trials_raw_diam_stims_reordered_entire=trials_raw_diam_stims_reordered_entire;

new_velocity_data.trials_raw_diam_difference_stimulation=trials_raw_diam_difference_stimulation;
new_velocity_data.trials_raw_diam_diff_relative_stimulation=trials_raw_diam_diff_relative_stimulation;
new_velocity_data.av_raw_diam_diff_relative_stimulation=av_raw_diam_diff_relative_stimulation;
    
new_velocity_data.av_prct_diam_stims_reordered_onset=av_prct_diam_stims_reordered_onset;
new_velocity_data.trials_prct_diam_stims_reordered_onset=trials_prct_diam_stims_reordered_onset;
new_velocity_data.av_prct_diam_stims_reordered_offset=av_prct_diam_stims_reordered_offset;
new_velocity_data.trials_prct_diam_stims_reordered_offset=trials_prct_diam_stims_reordered_offset;
new_velocity_data.av_prct_diam_stims_reordered_middle=av_prct_diam_stims_reordered_middle;
new_velocity_data.trials_prct_diam_stims_reordered_middle=trials_prct_diam_stims_reordered_middle;
new_velocity_data.av_prct_diam_stims_reordered_entire=av_prct_diam_stims_reordered_entire;
new_velocity_data.trials_prct_diam_stims_reordered_entire=trials_prct_diam_stims_reordered_entire;

new_velocity_data.trials_prct_diam_difference_stimulation=trials_prct_diam_difference_stimulation;

save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_velocity_data.mat','new_velocity_data')
