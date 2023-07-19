%% plot differences delta pupil
%% do calculation of still and running trials
% do reordering of trials so it starts from 0 deg - 330; calculates number
% of running trials depends on stimulus
% calculates pupil diameter for running and non-running trials
% calculates correlation between velocity and pupil diameter
% KS 2018/05/04
% open plot_peak_eye_depends_stim to create trials_eye_stims



for iAn=1:size(diam,2)
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
%eye_tmp=diam_norm{iAn}';
stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
epochs{iAn}=expt2{iAn}.frames.epochs;
stims{iAn}=expt2{iAn}.frames.stims;
%eye_tmp=diam_norm{iAn}; % looks better on the plots
%eye_tmp=zscore(diam{iAn}); % looks better on the plots
eye_tmp=diam{iAn}; % looks better on the plots
%eye_tmp=diam_prct_norm{iAn};
%eye_tmp=diam_stand{iAn}; % looks better on the plots
%eye_tmp=diam_zscore{iAn};
epochs_temp=expt2{iAn}.frames.epochs;

[av_eye_epochs{iAn} trials_eye_epochs{iAn}]=tcEpochAverage2(eye_tmp,epochs_temp);% time missing 

stims_temp=expt2{iAn}.frames.stims;
[av_eye_stims{iAn} trials_eye_stims{iAn}]=tcEpochAverage2(eye_tmp,stims_temp);% time missing 

clear eye_tmp
clear epochs_temp
clear stims_temp
end
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
%stim_eye_in=trials_eye_epochs{iAn};
stim_eye_in=trials_eye_stims{iAn};

frameRate=expt2{iAn}.frameRate;
time_epochs{iAn}=(0:(-1+size(stim_eye_in,1)))/frameRate; % time stamps

for iStim=1:expt2{iAn}.info.nStim-1;

[loc_trial loc_stim_tmp]=find(loctrial{iAn}(:,iStim)==1);
[still_trial still_stim_tmp]=find(stilltrial{iAn}(:,iStim)==1);
%[nondef_trial nondef_stim_tmp]=find(nondefined{iAn}(:,iStim)==1);
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
figure,

ax1_delta=axes('position',[0.1, 0.2, 0.27,0.16])
shadedErrorBar(1:size(eye_trace_alltrials,2),squeeze(nanmean(eye_trace_alltrials,1)),...
    nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'k');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
shadedErrorBar(1:size(eye_trace_still,2),squeeze(nanmean(eye_trace_still,1)),...
    nanstd(eye_trace_still,[],1)./sqrt(size(eye_trace_still,1)),'k');

ylabel('Diff Pupil area all trials(mm^2)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*4);
ylim([-0.02 0.1])

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'figure1_deltapupil_all_stationary_40expt.pdf']);

%% MEAN + SEM
figure
subplot(221)

errorbar(1:size(eye_trace_still,2),squeeze(nanmean(eye_trace_still,1)),...
      nanstd(eye_trace_still,[],1)./sqrt(size(eye_trace_still,1)),'k');

hold on
axis tight
ylim([-0.05 0.15]);
xlim([-0.03 13]);
ylabel('delta pupil size Stationary');
xlabel('Direction (deg)');
%title(sprintf('trial=%d',itrial));
set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(gca,'ticklength')*4);  

hold on
errorbar(1:size(eye_trace_locomotion,2),squeeze(nanmean(eye_trace_locomotion,1)),...
      nanstd(eye_trace_locomotion,[],1)./sqrt(size(eye_trace_locomotion,1)),'r');


subplot(222)

errorbar(1:size(eye_trace_still,2),squeeze(nanmean(eye_trace_still,1)),...
      nanstd(eye_trace_still,[],1)./sqrt(size(eye_trace_still,1)),'k');

hold on
axis tight
ylim([-0.05 0.15]);
xlim([-0.03 13]);
ylabel('delta pupil size Stationary');
xlabel('Direction (deg)');
%title(sprintf('trial=%d',itrial));
set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(gca,'ticklength')*4);  
%%

%% LOPCOMOTION VS STILL TRIALS
figure,
clf
clear tt
clear list_expt
for iAn=1:40; tt(iAn)= sum(eye_trace_still(iAn,:)); end
list_expt=find(tt>0);

ax1_delta=axes('position',[0.1, 0.2, 0.27,0.16])

errorbar(1:size(eye_trace_alltrials,2),squeeze(nanmean(eye_trace_alltrials,1)),...
      nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'b');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
errorbar(1:size(eye_trace_still(list_expt,:),2),squeeze(nanmean(eye_trace_still(list_expt,:),1)),...
      nanstd(eye_trace_still(list_expt,:),[],1)./sqrt(size(eye_trace_still(list_expt,:),1)),'k');

ylabel('Diff Pupil area all trials(mm^2)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*3);
ylim([-0.02 0.1])

clear tt
clear list_expt
for iAn=1:40; tt(iAn)= sum(eye_trace_locomotion(iAn,:)); end
list_expt=find(tt>0);

%ax1_delta=axes('position',[0.1, 0.5, 0.27,0.16])
errorbar(1:size(eye_trace_locomotion(list_expt,:),2),squeeze(nanmean(eye_trace_locomotion(list_expt,:),1)),...
      nanstd(eye_trace_locomotion(list_expt,:),[],1)./sqrt(size(eye_trace_locomotion(list_expt,:),1)),'r');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on

ylabel('Diff Pupil area all trials(mm^2)');
ylim([-0.02 0.1])

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'figure1_deltapupil_locomotion_all_stationary_40expt.pdf']);
