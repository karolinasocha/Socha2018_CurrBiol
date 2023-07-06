clear all

load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;

%%
clear av_eye_epochs
clear trials_eye_epochs
clear trials_eye_stims
clear av_eye_stims

%%
for iAn=1:size(diameter,2)
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
%eye_tmp=diam_norm{iAn}';
stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
epochs{iAn}=expt2{iAn}.frames.epochs;
stims{iAn}=expt2{iAn}.frames.stims;
%eye_tmp=diam_norm{iAn}; % looks better on the plots
%eye_tmp=zscore(diam{iAn}); % looks better on the plots
eye_tmp=diameter{iAn}; % looks better on the plots
%eye_tmp=diam_prct_norm{iAn};
%eye_tmp=diam_stand{iAn}; % looks better on the plots
%eye_tmp=diam_zscore{iAn};

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

%%

addpath 'G:\mousebox\code\mouselab\users\karolina\Attentional_Engagement_Locomotion'
clear time_epochs
clear resp_nondef_eye_trials
clear resp_still_eye_trials
clear resp_loc_eye_trials

    for iAn=1:size(velocity,2);

vtime=1:length(velocity{iAn});
%vval=velocity{iAn}(:,2);
vval=velocity{iAn};
locthresh=1;
stillthresh=1;
perthresh=0.95;
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

[loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials(1:length(vval),vval,...
    tstart,tstop,locthresh,stillthresh,perthresh);
clear stim_eye_in
nondefined{iAn}= ~(stilltrial{iAn}(:,:)+loctrial{iAn}(:,:));

%stim_eye_in=trials_eye_epochs{iAn};
stim_eye_in=trials_eye_stims{iAn};

frameRate=expt2{iAn}.frameRate;
time_epochs{iAn}=(0:(-1+size(stim_eye_in,1)))/frameRate; % time stamps

for iStim=1:expt2{iAn}.info.nStim-1;

[loc_trial loc_stim_tmp]=find(loctrial{iAn}(:,iStim)==1);
[still_trial still_stim_tmp]=find(stilltrial{iAn}(:,iStim)==1);
[nondef_trial nondef_stim_tmp]=find(nondefined{iAn}(:,iStim)==1);
% select pupil diameter eye after standarization without reordering
% stimulus
resp_loc_eye_trials{iAn}(:,iStim)=nan(size(time_epochs{iAn}))';
try
resp_loc_eye_trials{iAn}(:,iStim)=squeeze(nanmean(stim_eye_in(:,:,loc_trial',iStim),3));
end

resp_still_eye_trials{iAn}(:,iStim)=nan(size(time_epochs{iAn}))';
try
resp_still_eye_trials{iAn}(:,iStim)=squeeze(nanmean(stim_eye_in(:,:,still_trial',iStim),3));
end

resp_nondef_eye_trials{iAn}(:,iStim)=nan(size(time_epochs{iAn}))';
try
resp_nondef_eye_trials{iAn}(:,iStim)=squeeze(nanmean(stim_eye_in(:,:,nondef_trial',iStim),3));
end


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
    end

%%
j=0;

eye_trace_still=nan([size(newdir,1) 12]); 
eye_trace_alltrials=nan([size(newdir,1) 12]);
for iAn=1:size(newdir,1)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2{iAn}.frameRate);
        epochs_length=length(resp_still_eye_trials{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

eye_trace_still(iAn,:)=squeeze(nanmean(resp_still_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_still_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     

end

%%
%%
j=0;

eye_trace_locomotion=nan([size(newdir,1) 12]); 
eye_trace_alltrials=nan([size(newdir,1) 12]);
for iAn=1:size(newdir,1)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2{iAn}.frameRate);
        epochs_length=length(resp_loc_eye_trials{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

eye_trace_locomotion(iAn,:)=squeeze(nanmean(resp_loc_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_loc_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     

end
%% MIXED

%%
j=0;

eye_trace_nondef=nan([size(newdir,1) 12]); 
eye_trace_alltrials=nan([size(newdir,1) 12]);
for iAn=1:size(newdir,1)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2{iAn}.frameRate);
        epochs_length=length(resp_nondef_eye_trials{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

eye_trace_nondef(iAn,:)=squeeze(nanmean(resp_nondef_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_nondef_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     

end

%% LOPCOMOTION VS STILL TRIALS
figure,
clear tt
clear list_expt
for iAn=1:40; tt(iAn)= sum(eye_trace_still(iAn,:)); end
list_expt=find(tt>0);

ax1_delta=axes('position',[0.1, 0.2, 0.27,0.16])
shadedErrorBar(1:size(eye_trace_alltrials,2),squeeze(nanmean(eye_trace_alltrials,1)),...
    nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'k');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
shadedErrorBar(1:size(eye_trace_still(list_expt,:),2),squeeze(nanmean(eye_trace_still(list_expt,:),1)),...
    nanstd(eye_trace_still(list_expt,:),[],1)./sqrt(size(eye_trace_still(list_expt,:),1)),'b');

ylabel('Diff Pupil area all trials(mm^2)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*4);
ylim([-0.02 0.1])


clear tt
clear list_expt
for iAn=1:40; tt(iAn)= sum(eye_trace_locomotion(iAn,:)); end
list_expt=find(tt>0);

ax1_delta=axes('position',[0.1, 0.5, 0.27,0.16])
shadedErrorBar(1:size(eye_trace_alltrials,2),squeeze(nanmean(eye_trace_alltrials,1)),...
    nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'k');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
shadedErrorBar(1:size(eye_trace_locomotion(list_expt,:),2),squeeze(nanmean(eye_trace_locomotion(list_expt,:),1)),...
    nanstd(eye_trace_locomotion(list_expt,:),[],1)./sqrt(size(eye_trace_locomotion(list_expt,:),1)),'r');

ylabel('Diff Pupil area all trials(mm^2)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*4);
ylim([-0.02 0.1])


%%
clear tt
clear list_expt
for iAn=1:40; tt(iAn)= sum(eye_trace_nondef(iAn,:)); end
list_expt=find(tt>0);

ax1_delta=axes('position',[0.1, 0.75, 0.27,0.16])
shadedErrorBar(1:size(eye_trace_alltrials,2),squeeze(nanmean(eye_trace_alltrials,1)),...
    nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'k');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
shadedErrorBar(1:size(eye_trace_nondef(list_expt,:),2),squeeze(nanmean(eye_trace_nondef(list_expt,:),1)),...
    nanstd(eye_trace_nondef(list_expt,:),[],1)./sqrt(size(eye_trace_nondef(list_expt,:),1)),'--r');

ylabel('Diff Pupil area all trials(mm^2)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*4);
ylim([-0.02 0.1])

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'figure1_deltapupil_stationary_locomotion_all_40expt.pdf']);

%% 
clear eye_array_alltrials
all_array_trials=nan(40,15,12);
size(all_array_trials)

for iAn=1:40
        clear correct_value
        clear order
        frameRate=round(expt2{iAn}.frameRate);
        epochs_length=length(resp_nondef_eye_trials{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end
        
        clear tmp_data
        tmp_data=squeeze(trials_eye_stims{iAn});
        eye_array_alltrials{iAn}(:,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,:,order)))-...
        squeeze(nanmean(tmp_data(frameRate:2*frameRate,:,order)));     
        all_array_trials(iAn,1:size(tmp_data,2),:)=eye_array_alltrials{iAn};
end
      
for itrial=1:15
eye_trial(itrial,:)=squeeze(nanmean(all_array_trials(:,itrial,:),1));
end

figure
colors_trials=hsv(10);

for itrial=1:10
plot(eye_trial(itrial,:),'color',colors_trials(itrial,:));
hold on
end

size(all_array_trials)
%% PLOT ALL TRIALS
figure
colors_trials=hsv(10);
ax=[];
for itrial=1:10
%ax(itrial)=subplot(5,2,itrial);
subplot(222);
hold on
clear temporary_data
temporary_data=squeeze(all_array_trials(:,itrial,:));
h1=shadedErrorBar(1:size(temporary_data,2),squeeze(nanmean(temporary_data,1)),...
    nanstd(temporary_data(:,:),[],1)./sqrt(size(temporary_data,1)),'k',0.5);
%set(gca,'color',colors_trials(itrial,:));
h1.patch.FaceColor=colors_trials(itrial,:);
hold on
axis tight


% set(ax(itrial),'xtick',[1:3:12],'xticklabel',[0:90:330],...
%     'tickdir','out','box','off','layer','top','color','none',...
%     'fontsize',12,'ticklength',get(ax(itrial),'ticklength')*4);
ylim([-0.02 0.12]);
%title(sprintf('trial=%d',itrial));

end

set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(gca,'ticklength')*4);

%%
   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_diameter_trials\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'mergedplot_deltapupil_trials_average_40expt.pdf']);

%% PLOT ALL TRIALS
figure
colors_trials=hsv(10);
ax=[];
err=[];
for itrial=1:10
%ax(itrial)=subplot(5,2,itrial);
subplot(222);
hold on
clear temporary_data
temporary_data=squeeze(all_array_trials(:,itrial,:));
err(itrial)=errorbar(1:size(temporary_data,2),squeeze(nanmean(temporary_data,1)),...
    nanstd(temporary_data(:,:),[],1)./sqrt(size(temporary_data,1)));
set(err(itrial),'color',[0.5 0.5 0.5]);

if itrial==2;
set(err(itrial),'color',[1 0 0],'linewidth',1.5);
end

if itrial==10;
set(err(itrial),'color',[0 0 0],'linewidth',1.5);
end
%colors_trials(itrial,:)
hold on
axis tight

% set(ax(itrial),'xtick',[1:3:12],'xticklabel',[0:90:330],...
%     'tickdir','out','box','off','layer','top','color','none',...
%     'fontsize',12,'ticklength',get(ax(itrial),'ticklength')*4);
ylim([-0.05 0.15]);
xlim([-0.03 13]);
ylabel('delta pupil size');
xlabel('Direction (deg)');
%title(sprintf('trial=%d',itrial));

end

set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(gca,'ticklength')*4);
  
set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_diameter_trials\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'comparison_deltapupil_trials_average_40expt.pdf']);
%%
figure
colors_trials=hsv(10);
ax=[];
for itrial=1:10
ax(itrial)=subplot(5,2,itrial);
hold on
clear temporary_data
temporary_data=squeeze(all_array_trials(:,itrial,:));
h1=shadedErrorBar(1:size(temporary_data,2),squeeze(nanmean(temporary_data,1)),...
    nanstd(temporary_data(:,:),[],1)./sqrt(size(temporary_data,1)),'k',0.5);
set(gca,'color',colors_trials(itrial,:));
h1.patch.FaceColor=colors_trials(itrial,:);
hold on
axis tight

set(ax(itrial),'xtick',[1:3:12],'xticklabel',[0:90:330],...
     'tickdir','out','box','off','layer','top','color','none',...
     'fontsize',12,'ticklength',get(ax(itrial),'ticklength')*4);
ylim([-0.02 0.12]);
title(sprintf('trial=%d',itrial));

end

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_diameter_trials\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'singleplots_deltapupil_trials_average_40expt.pdf']);
%%

