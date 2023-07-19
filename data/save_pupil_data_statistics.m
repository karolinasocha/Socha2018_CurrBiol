% create zscored pupil_data for statistical tests

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
clear av_eye_epochs
clear trials_eye_epochs
clear trials_eye_stims
clear av_eye_stims

%% creates epochs for evry stimulus condition and trial; 
% the length depends on the framerate of acquisition so it matches imaging
% data

for iAn=1:size(diameter,2)
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
epochs{iAn}=expt2{iAn}.frames.epochs;
stims{iAn}=expt2{iAn}.frames.stims;
%eye_tmp=diam{iAn}; % looks better on the plots
diam_norm_median{iAn}=(diam{iAn}-nanmedian(diam{iAn}))./nanmedian(diam{iAn});
eye_tmp=diam_norm_median{iAn};

epochs_temp=expt2{iAn}.frames.epochs;
try
[av_eye_epochs{iAn} trials_eye_epochs{iAn}]=tcEpochAverage2(eye_tmp,epochs_temp);% time missing 
end
try
[av_eye_epochs{iAn} trials_eye_epochs{iAn}]=tcEpochAverage2(eye_tmp',epochs_temp);% time missing 
end
stims_temp=expt2{iAn}.frames.stims;
try
[av_eye_stims{iAn} trials_eye_stims{iAn}]=tcEpochAverage2(eye_tmp,stims_temp);% time missing 
end
try
[av_eye_stims{iAn} trials_eye_stims{iAn}]=tcEpochAverage2(eye_tmp',stims_temp);% time missing 
end
clear eye_tmp
clear epochs_temp
clear stims_temp
end

%%
%% creates eye_trace_alltrials
% average zscored pupil during visual stimulation for individual trials

j=0;
eye_trace_alltrials=nan([size(newdir,1) 12]);

for iAn=1:size(newdir,1)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2{iAn}.frameRate);
        epochs_length=length(trials_eye_stims{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

clear tmp_data

%tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
for itrial=1:size(trials_eye_stims{iAn},3);
    tmp_data=squeeze(trials_eye_stims{iAn}(:,:,itrial,:));
%eye_trace_individualtrials{iAn}(itrial,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
%squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     
%eye_trace_individualtrials{iAn}(itrial,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)));     
eye_trace_individualtrials{iAn}(itrial,:)=squeeze(nanmean(tmp_data(:,order)));     
eye_trace_individualtrials_control{iAn}(itrial)=squeeze(nanmean(tmp_data(:,13)));

end
end

%% total size
for iAn=1:length(eye_trace_individualtrials)
    num_trials(iAn)=size(eye_trace_individualtrials{iAn},1);
end

%%
pupil_nan_array = NaN(sum(num_trials), 16)

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=eye_trace_individualtrials{iAn};
    tmp_control=eye_trace_individualtrials_control{iAn};
    
    for itrial=1:size(tmp,1)
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);
        pupil_nan_array(t,16)=tmp_control(itrial);
    end
end

%%
pupil_data.directions=correct_value
pupil_data.animal_id=animal_id;
pupil_data.sessions_id=session_id;

pupil_data.zscored_data=pupil_nan_array;

save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\pupil_data.mat','pupil_data')
%%
nconditions=13

data_all={}; % cells containing the parameters of inidividual areas
for icond=1:nconditions
    data_all{icond}=pupil_data.zscored_data(:,icond+3);
end
%%
threshold=0.05; % default significance level 95%
cond_name={'0','30','60','90','120','150','180','210','240','270','300','330','control'}
options.nruns = 1000; % # random sampling
options.nanimal = 5;  % # mice per sampling
options.npopsize = 100; % # neurons per mouse
options.isplot = 0; % plot intermediate steps or not

figure(2)
clf
plot_significance_bootstrpKS(data_all,[],options); % by default, from light to dark blue, p>0.05, <0.05, <0.01, <0.001
set(gca,'xtick',[1:nconditions],'xticklabel',cond_name);
set(gca,'ytick',[1:nconditions],'yticklabel',cond_name);
xtickangle(90);

