
%% this is script to generate Figure3F Stationary trials comparison
% edited 02-August-2023 by KSocha
% it has weighted probability to select session

%% plot density plots
%% load data
% in relation to pupil diameter
clear all
%% plot pupil response median vs fraction of running
% G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\recalculation_pupil_median\figures\area_dynamics_behavior_condidtions
% pupil_alltrials_alltrials_fractionrunning_average_pupil.pdf

clear all

%load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
newdir=new_pupil_data.newdir;
stimulus=new_pupil_data.stimulus;
velocity=new_pupil_data.velocity;
expt=new_pupil_data.expt;
expt2=new_pupil_data.expt2;

diameter_data=new_pupil_data.diameter
diameter_data{10}=diameter_data{10}';


%% it needs delta pupil
trials_raw_diam_diff_relative_stimulation=new_pupil_data.trials_raw_diam_diff_relative_stimulation;
raw_relative_diam_delta=new_pupil_data.raw_relative_diam_delta;
%% VELOCITY
array_data_velo=[];
matrix_data_velo=[];
velocity_data=velocity;

for iAn=1:size(newdir,1);
    clear av
    clear tmp1
    stims{iAn}=expt2{iAn}.frames.stims;  
    [av tmp1]=tcEpochAverage2(velocity_data{iAn},stims{iAn});
    trials_velocity_stims{iAn}=tmp1;
clear tmp_file
clear tmp_file
for iStim=1:size(tmp1,4)
    dtime = linspace(0,1,length(tmp1));
    p2time =  linspace(0,1,159);
tmp_file(:,:,iStim)=interp1(dtime,squeeze(tmp1(:,:,:,iStim)),p2time,'linear');

end
array_data_velo{iAn}=tmp_file;
matrix_data_velo{iAn}=tmp_file(:,:);
end
%%
%%
for iAn=1:length(newdir)
    newdir_animals_id{iAn}=newdir{iAn}(8:12);
end
animal_id_list=unique(newdir_animals_id);
animal_id=nan(length(newdir),1);
sessions_id=nan(length(newdir),1);

for iAn=1:length(animal_id_list)
    indexes=find(strcmp(newdir_animals_id, animal_id_list{iAn}));
    animal_id(indexes)=iAn;
    % session
    sessions_id(indexes)=1:length(indexes)
end

% session

%%
velocity_data=velocity;
for iAn=1:size(velocity_data,2)
        epochs{iAn}=expt2{iAn}.frames.epochs;
        stims{iAn}=expt2{iAn}.frames.stims;
        tc_boutons{iAn}=diameter_data{iAn};
        [tc_response_epochs{iAn} tc_trial_response_epochs{iAn}]=tcEpochAverage2(tc_boutons{iAn},epochs{iAn});
        [tc_response_stims{iAn} tc_trial_response_stims{iAn}]=tcEpochAverage2(tc_boutons{iAn},stims{iAn});    
        [tc_velocity_stims{iAn} tc_trial_velocity_stims{iAn}]=tcEpochAverage2(velocity{iAn},stims{iAn});
        [tc_velocity_epochs{iAn} tc_trial_velocity_epochs{iAn}]=tcEpochAverage2(velocity{iAn},epochs{iAn});
end
 
%%
clear response_norm_baseline
for iAn=1:size(tc_boutons,2)
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
%eye_tmp=diam{iAn};

stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
epochs{iAn}=expt2{iAn}.frames.epochs;
stims{iAn}=expt2{iAn}.frames.stims;
epochs_temp=expt2{iAn}.frames.epochs;
%[av_eye_epochs{iAn} trials_eye_epochs{iAn}]=tcEpochAverage2(eye_tmp,epochs_temp);% time missing 

stims_temp=expt2{iAn}.frames.stims;
%[av_eye_stims{iAn} trials_eye_stims{iAn}]=tcEpochAverage2(eye_tmp,stims_temp);% time missing 

[av_velo{iAn} trials_velocity_stims{iAn}]=tcEpochAverage2(velocity{iAn},stims_temp);
clear eye_tmp
clear epochs_temp
clear stims_temp
clear stimulus_resp_tmp

end

%%
clear resp_loc_velo_trials
clear resp_still_velo_trials
clear resp_nondef_velo_trials
clear resp_still_eye_trials 
clear resp_nondef_eye_trials
clear resp_loc_eye_trials
clear loctrial
clear stilltrial

%
addpath 'G:\mousebox\code\mouselab\users\karolina\Attentional_Engagement_Locomotion'

    for iAn=1:size(velocity,2);

vtime=1:length(velocity{iAn});
vval=smoothdata(velocity{iAn},'sgolay'); % used to be smooth
clear tstart
clear tstop
% locthresh=1;
% stillthresh=0.5;
% perthresh_stationary=0.95;
% perthresh=0.05;

locthresh=1;
stillthresh=1;
perthresh_stationary=0.95;
perthresh=0.05;

tstart = nan(size(expt2{iAn}.frames.stims));
k=0;
for i=1:expt2{iAn}.info.nTrials; 
    for j=1:expt2{iAn}.info.nStim;
        k=k+1;
        tstart(i,j)=expt2{iAn}.frames.stims{i,j}(1); %expt2{iAn}.frames.stims{i,j}(1);
        %tstart_vec(k)=expt2{iAn}.frames.stims{i,j}(1)
        tstop(i,j)=expt2{iAn}.frames.stims{i,j}(end);
        %tstop_vec(k)=expt2{iAn}.frames.stims{i,j}(end)
    end; 
end;

% [loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials(1:length(vval),vval,...
%     tstart,tstop,locthresh,stillthresh,perthresh);

  [loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials_val_stationary(1:length(vval),vval,...
      tstart,tstop,locthresh,stillthresh,perthresh,perthresh_stationary);

[loc_trial loc_stim_tmp]=find(loctrial{iAn}(:,1:end-1)==1);
[still_trial still_stim_tmp]=find(stilltrial{iAn}(:,1:end-1)==1);
nondefined{iAn}= ~(stilltrial{iAn}(:,:)+loctrial{iAn}(:,:));
[nondef_trial nondef_stim_tmp]=find(nondefined{iAn}(:,1:end-1)==1);
% select pupil diameter eye after standarization without reordering
% stimulus

tc_in=tc_trial_response_stims{iAn};

resp_loc_resp_trials{iAn}=tc_in(:,:,loc_trial',loc_stim_tmp');
resp_still_resp_trials{iAn}=tc_in(:,:,still_trial',still_stim_tmp');
resp_nondef_resp_trials{iAn}=tc_in(:,:,nondef_trial',nondef_stim_tmp');

%
stim_deg_loc{iAn}=stimulus{iAn}(loc_stim_tmp);
stim_deg_still{iAn}=stimulus{iAn}(still_stim_tmp);
stim_deg_nondef{iAn}=stimulus{iAn}(nondef_stim_tmp);

clear vtime
clear vval
clear tstart
clear tstop
clear tc_loc_tmp
clear tc_still_tmp
clear tc_nondef_tmp
clear loc_stim_tmp
clear still_stim_tmp
clear nondef_stim_tmp
clear loc_trial
clear still_trial
clear nondef_trial
    end
%%
for iAn=1:length(loctrial)
    ntrials_locomotion(iAn,:)=sum(loctrial{iAn});
    ntrials_still(iAn,:)=sum(stilltrial{iAn});
    tmp_still=ntrials_still(iAn,1:12);
    tmp_loc=ntrials_locomotion(iAn,1:12)
    if sum(length(find(tmp_still>2)))==12
        session_still(iAn)=1
    else
        session_still(iAn)=0
    end
    
    if sum(length(find(tmp_loc>2)))==12
        
        session_loc(iAn)=1
    else
        session_loc(iAn)=0
    end
                 
end
%%
