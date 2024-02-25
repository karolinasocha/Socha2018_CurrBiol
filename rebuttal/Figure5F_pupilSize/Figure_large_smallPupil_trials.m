% % from Sylvia Schroeder Methods
% Two different tuning curves were fitted for responses during small and large pupil trials. To categorize each trial, the trace of pupil size
% was interpolated to match the sampling times of the neural data (spikes were binned into 5 ms bins) 
% and then thresholded at the median pupil size measured during the experiment. 
% If pupil size was below that threshold for longer than half of the trial, the trial was
% categorized as small pupil trial, and otherwise as large pupil trial.

clear all
% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
% new_pupil_data.trials_raw_diam_stims_reordered_onset

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
diam=diameter;
session_id= data_behav_40expt.session_id;
animal_id= data_behav_40expt.animal_id;

%%
%% creates epochs and re-order stimulus so it is consistent for all experiments 

for iAn=1:size(diameter,2)


diam_raw=diam{iAn};
pupil_median(iAn)=nanmedian(diam{iAn})
all_epochs=expt2{iAn}.frames.epochs;
stim_epochs=expt2{iAn}.frames.stims;

% raw data
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_raw;
[av_raw_diam_epochs{iAn},trials_raw_diam_epochs{iAn}, av_raw_diam_stims{iAn}, trials_raw_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_raw_diam_epochs_reordered{iAn}, trials_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_raw_diam_epochs{iAn},trials_raw_diam_stims{iAn});
[av_raw_diam_epochs_reordered{iAn}, av_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_raw_diam_epochs{iAn},av_raw_diam_stims{iAn});

for iTrial=1:size(trials_raw_diam_stims_reordered{iAn},3)
for iStim=1:size(trials_raw_diam_stims_reordered{iAn},4)
    fraction_smallPupil{iAn}(iTrial,iStim)=sum(trials_raw_diam_stims_reordered{iAn}(:,:,iTrial,iStim)<pupil_median(iAn))/size(trials_raw_diam_stims_reordered{iAn},1)
    if fraction_smallPupil{iAn}(iTrial,iStim)>0.5
        trial_smallPupil{iAn}(iTrial,iStim)=1
    else
        trial_smallPupil{iAn}(iTrial,iStim)=0
    end
    fraction_largePupil{iAn}(iTrial,iStim)=sum(trials_raw_diam_stims_reordered{iAn}(:,:,iTrial,iStim)>pupil_median(iAn))/size(trials_raw_diam_stims_reordered{iAn},1)
    if fraction_largePupil{iAn}(iTrial,iStim)>=0.5
        trial_largePupil{iAn}(iTrial,iStim)=1
    else
        trial_largePupil{iAn}(iTrial,iStim)=0
    end
end
end

end
%% proportion of number of trials with large pupil vs small pupil

ntrials_largePupil=sum(cell2mat(trial_largePupil'))/length(cell2mat(trial_largePupil'))
ntrials_smallPupil=sum(cell2mat(trial_smallPupil'))/length(cell2mat(trial_largePupil'))
%%
figure, 
plot(ntrials_smallPupil,'k','linewidth',2)
hold on
plot(ntrials_largePupil,'r','linewidth',2)
ylabel('Fraction trials','fontsize',18)
xlabel('Stimulus (deg)','fontsize',18)
legend('Small Pupil', 'Large Pupil', 'Location', 'best'); % Add a legend
set(gca, 'TickDir', 'out'); % Set ticks to point outwards
set(gca, 'XTick', [1:1:13]); % Set custom x-tick positions
set(gca, 'XTickLabel', [0,30,60,90,120,150,180,210,240,270,300,330]); % Set custom x-tick labels
box off


set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\Figure5F_pupilSize\']; 
filename=sprintf(['Figure_same_large_pupil_fractionTrials',newdir_awake(1:18),'.pdf'])

print(gcf,'-dpdf',[filepathanalysis, filename]);


%%
for iAn=1:length(trial_smallPupil)
ntrials_small_per_animal(iAn,:)=sum(trial_smallPupil{iAn})
ntrials_large_per_animal(iAn,:)=sum(trial_largePupil{iAn})
end

%%

ntrials_small_per_animal


%% creates epochs and re-order stimulus so it is consistent for all experiments 

for iAn=1:size(diameter,2)

pupil_median(iAn)=nanmedian(diam{iAn})
trials_epochs=new_pupil_data.av_raw_diam_stims_reordered_onset{iAn}

end
%%
trials_epochs=new_pupil_data.av_raw_diam_stims_reordered_onset{iAn}
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
diam=diameter;
session_id= data_behav_40expt.session_id;
animal_id= data_behav_40expt.animal_id;
%%
figure, plot(diam{iAn})
%%
iAn=1
diam_norm_median=(diam{iAn}-nanmedian(diam{iAn}))./nanmedian(diam{iAn});
diam_example=diam{iAn}
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
