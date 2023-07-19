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

end
end

%% averages zscored pupil size across trials
for iAn=1:40
eye_trace_individualtrials_av{iAn}=nanmean(eye_trace_individualtrials{iAn},1)
end
%% creates matrix all sessions and mice

nSess=size(newdir,1)
[nT nS]=size(eye_trace_individualtrials{1})
session_array=nan(nSess,nT,nS)

for nSess=1:size(newdir,1);
[nT nS]=size(eye_trace_individualtrials{nSess})

session_array(nSess,1:nT,1:nS)=eye_trace_individualtrials{nSess};
session_array_aver(nSess,1:nS)=eye_trace_individualtrials_av{nSess};

end

%% plot same animals SEM

addpath 'G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals'
pattern_string='KS'
[animals_array]=find_sameanimals(newdir,pattern_string)
[x y]=find(animals_array==iAn)
%animals_array(iAn,~isnan(animals_array(iAn,:)))

for iAn=1:40
    [x y]=find(animals_array==iAn);
    x_animal(iAn)=x
    y_session(iAn)=y
end

%% arbitrary distinguish which directions is nasal vs temporal
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];
clear eye_trace_nasal_sameanimal
clear eye_trace_temporal_sameanimal

for n_animals=1:size(animals_array,1)
eye_trace_alltrials_sameanimal{n_animals}=...
    session_array(find(x_animal(:)==n_animals),:,:);

eye_trace_alltrials_sameanimal_aver{n_animals}=...
    session_array_aver(find(x_animal(:)==n_animals),:);

eye_trace_nasal_sameanimal{n_animals}=...
    session_array_aver(find(x_animal(:)==n_animals),order_nasal);

eye_trace_temporal_sameanimal{n_animals}=...
    session_array_aver(find(x_animal(:)==n_animals),order_temporal);
end
%%
clear trials_animals_temporal
clear trials_animals_nasal
for iAn=1:13
trace_animals(iAn,:)=squeeze(nanmean(eye_trace_alltrials_sameanimal_aver{iAn},1));
trace_animals_nasal(iAn,:)=squeeze(nanmean(eye_trace_nasal_sameanimal{iAn},1));
trace_animals_temporal(iAn,:)=squeeze(nanmean(eye_trace_temporal_sameanimal{iAn},1));

trials_animals_nasal{iAn}=eye_trace_alltrials_sameanimal{iAn}(:,:,order_nasal);
trials_animals_temporal{iAn}=eye_trace_alltrials_sameanimal{iAn}(:,:,order_temporal);

end
%% 
for iAn=1:13
    clear tmp1
    clear tmp2
    tmp1=trials_animals_temporal{iAn};
    tmp1(isnan(tmp1))=[];
    
    tmp2=trials_animals_nasal{iAn};
    tmp2(isnan(tmp2))=[];
    
    array_animals_temporal_trials{iAn}=tmp1;
    array_animals_nasal_trials{iAn}=tmp2;

end
%% kstest for all individual trials
for iAn=1:13
[h0_kstest_trials(iAn), pval_kstest_trials(iAn)]=kstest2(array_animals_temporal_trials{iAn}(:), array_animals_nasal_trials{iAn}(:));

end

%%
figure, 
imagesc(trace_animals)
% test statistics for each animal
clear p
clear h
for iAn=1:13
    clear tmp_nasal
    clear tmp_temporal
tmp_temporal=trace_animals_temporal(iAn,:);
tmp_nasal=trace_animals_nasal(iAn,:);

[h(iAn), p(iAn)]=ttest2(tmp_nasal, tmp_temporal);

%[h(iAn), p(iAn)]=kstest2(tmp_nasal, tmp_temporal);

end
colormap(bluered)
pval_kstest=p
h0_kstest=h
%% average mouse heat map
fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');

axes(ax_heatmap_allmice_plot)
imagesc(trace_animals)
axis tight
axis square

cmap=colormap(redblue)
set(gca,'CLim',[-0.4 0.4])

cb=colorbar;
cb.Position = [.27,.082,.015,.05] 
caxis([-0.2 0.3]);

clear average_trials_animal
for iAn=1:13
average_trials_animal(iAn,:,:)=squeeze(nanmean(eye_trace_alltrials_sameanimal{iAn},1));
end

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

for iAn=1:13
   av_nasal_trials(iAn,:,:)= average_trials_animal(iAn,:,order_nasal);
   av_temporal_trials(iAn,:,:)= average_trials_animal(iAn,:,order_temporal);   
end

for iAn=1:13;
    
    clear tmp_temporal
    clear tmp_nasal
    
tmp_temporal=av_temporal_trials(iAn,:);
tmp_nasal=av_nasal_trials(iAn,:);
%[h0(iAn), pval(iAn)]=ttest2(tmp_nasal, tmp_temporal);
[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(tmp_nasal(~isnan(tmp_nasal)), tmp_temporal(~isnan(tmp_temporal)), 'Tail','right');
%[pval(iAn),h0(iAn),stats{iAn}]=signrank(tmp_nasal(~isnan(tmp_nasal)), tmp_temporal(~isnan(tmp_temporal)), 'Tail','right');

%[h0(iAn), pval(iAn)]=kstest2(tmp_nasal, tmp_temporal);

end

[h0; pval]
for iAn=1:13
if pval(iAn)<0.05 & h0(iAn)==1
    str_name{iAn}=sprintf('* %d', iAn);    
else
    str_name{iAn}=sprintf('  %d', iAn);    
end
end
xtickstr={'0',[],[],'90',[],[],'180',[],[],'270',[],[]};
set(ax_heatmap_allmice_plot,'xtick',[1:1:12],'xticklabel',xtickstr,'ytick',[1:1:13],'yticklabel',str_name,...
     'tickdir','out','box','off','layer','top','color','none',...
     'fontsize',10,'ticklength',get(ax_heatmap_allmice_plot,'ticklength')*4);

set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\Figure1E\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_median_pupil_pvals_correctedcolorbar_2.pdf']);

%% heat map with kstest from individual trials (all trials)
fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');

axes(ax_heatmap_allmice_plot)
imagesc(trace_animals)
axis tight
axis square
cmap=colormap(redblue)
set(gca,'CLim',[-0.4 0.4])

cb=colorbar;
cb.Position = [.27,.082,.015,.05] 
caxis([-0.2 0.3]);

clear average_trials_animal
for iAn=1:13
average_trials_animal(iAn,:,:)=squeeze(nanmean(eye_trace_alltrials_sameanimal{iAn},1));
end

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

for iAn=1:13
   av_nasal_trials(iAn,:,:)= average_trials_animal(iAn,:,order_nasal);
   av_temporal_trials(iAn,:,:)= average_trials_animal(iAn,:,order_temporal);   
end


for iAn=1:13
if pval_kstest(iAn)<0.05 & h0_kstest(iAn)==1
    str_name{iAn}=sprintf('* %d', iAn);    
else
    str_name{iAn}=sprintf('  %d', iAn);    
end
end
xtickstr={'0',[],[],'90',[],[],'180',[],[],'270',[],[]};
set(ax_heatmap_allmice_plot,'xtick',[1:1:12],'xticklabel',xtickstr,'ytick',[1:1:13],'yticklabel',str_name,...
     'tickdir','out','box','off','layer','top','color','none',...
     'fontsize',10,'ticklength',get(ax_heatmap_allmice_plot,'ticklength')*4);

set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals\plot_heatmaps_pupil_median_pvals\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_median_pupil_pvals_correctedcolorbar_2.pdf']);



%%

for iAn=1:13;
tmp_temporal=trace_animals_temporal(iAn,:);
tmp_nasal=trace_animals_nasal(iAn,:);
[h0(iAn), pval(iAn)]=ttest2(tmp_nasal, tmp_temporal);
%[h0(iAn), pval(iAn)]=kstest2(tmp_nasal, tmp_temporal);

end

%% test statistics for each animal
clear p
clear h
for iAn=1:13
    clear tmp_nasal
    clear tmp_temporal
tmp_nasal=eye_trace_nasal_sameanimal{iAn}(:);
tmp_temporal=eye_trace_temporal_sameanimal{iAn}(:);

tmp_temporal(isnan(tmp_nasal(:)),:) = [];
tmp_nasal(isnan(tmp_nasal(:)),:) = [];
[h(iAn), p(iAn)]=ttest2(tmp_nasal, tmp_temporal);

[h(iAn), p(iAn)]=kstest2(tmp_nasal, tmp_temporal);

end

%%
clear av_animal
for iAn=1:13
av_animal(iAn,:)=nanmean(nanmean(eye_trace_alltrials_sameanimal{iAn},1),2);
end

clear h0
clear pval

for iAn=1:13
    clear tmp_nasal
    clear tmp_temporal
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];
tmp_nasal=av_animal(iAn,order_nasal);
tmp_temporal=av_animal(iAn,order_temporal);

[h0(iAn), pval(iAn)]=kstest2(tmp_nasal, tmp_temporal);

end



%% PLOT ONLY 0-30-330
clear all

load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;
%%
j=0;
for iAn=1:size(newdir,1)
   % SELECT EXPERIMENTS 
%if stimulus{iAn}(1)==0; 
    j=j+1;
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
%eye_tmp=diam_norm{iAn}';
newdir_sel{j}=newdir{iAn};
diameter_sel{j}=diameter{iAn};
expt_sel{j}=expt{iAn};
expt2_sel{j}=expt2{iAn};
velocity_sel{j}=velocity{iAn};
stimulus_sel{j}=[expt2{iAn}.prot.pars.ori];
epochs_sel{j}=expt2{iAn}.frames.epochs;
stims_sel{j}=expt2{iAn}.frames.stims;
%eye_tmp=diam_norm{iAn}; % looks better on the plots
%eye_tmp=zscore(diam{iAn}); % looks better on the plots
eye_tmp=diameter{iAn}; % looks better on the plots
%eye_tmp=diam_prct_norm{iAn};
%eye_tmp=diam_stand{iAn}; % looks better on the plots
%eye_tmp=diam_zscore{iAn};

epochs_temp=expt2{iAn}.frames.epochs;
try
[av_eye_epochs{j} trials_eye_epochs{j}]=tcEpochAverage2(eye_tmp',epochs_temp);% time missing 
end

try
[av_eye_epochs{j} trials_eye_epochs{j}]=tcEpochAverage2(eye_tmp,epochs_temp);% time missing 
end

stims_temp=expt2{iAn}.frames.stims;
try
[av_eye_stims{j} trials_eye_stims{j}]=tcEpochAverage2(eye_tmp',stims_temp);% time missing 
end

try
[av_eye_stims{j} trials_eye_stims{j}]=tcEpochAverage2(eye_tmp,stims_temp);% time missing 
end
clear eye_tmp
clear epochs_temp
clear stims_temp
%end

end
    
%%
clear av_eye_epochs
clear trials_eye_epochs
clear trials_eye_stims
clear av_eye_stims
%%
%%

%%

addpath 'G:\mousebox\code\mouselab\users\karolina\Attentional_Engagement_Locomotion'
clear time_epochs
clear resp_nondef_eye_trials
clear resp_still_eye_trials
clear resp_loc_eye_trials

    for iAn=1:size(velocity_sel,2);

        clear vval
        clear vtime
vtime=1:length(velocity_sel{iAn});
%vval=velocity{iAn}(:,2);
vval=velocity_sel{iAn};
locthresh=1;
stillthresh=1;
perthresh=0.95;
tstart = nan(size(expt2_sel{iAn}.frames.stims));
k=0;
for i=1:expt2_sel{iAn}.info.nTrials; 
    for j=1:expt2_sel{iAn}.info.nStim;
        k=k+1;
        tstart(i,j)=expt2_sel{iAn}.frames.stims{i,j}(1); %expt2{iAn}.frames.stims{i,j}(1);
        %tstart_vec(k)=expt2{iAn}.frames.stims{i,j}(1)
        tstop(i,j)=expt2_sel{iAn}.frames.stims{i,j}(end);
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

for iStim=1:expt2_sel{iAn}.info.nStim-1;

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

eye_trace_still=nan([size(newdir_sel,1) 12]); 
eye_trace_alltrials=nan([size(newdir_sel,1) 12]);
for iAn=1:size(newdir_sel,2)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2_sel{iAn}.frameRate);
        stims_value=stimulus_sel{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end
        
epochs_length=length(resp_still_eye_trials{iAn});
if epochs_length>0
eye_trace_still(iAn,:)=squeeze(nanmean(resp_still_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_still_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     
end

end

%%
%%
j=0;

eye_trace_locomotion=nan([size(newdir_sel,1) 12]); 
eye_trace_alltrials=nan([size(newdir_sel,1) 12]);
for iAn=1:size(newdir_sel,2)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2_sel{iAn}.frameRate);
        epochs_length=length(resp_loc_eye_trials{iAn});
        stims_value=stimulus_sel{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

        if epochs_length>0
eye_trace_locomotion(iAn,:)=squeeze(nanmean(resp_loc_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_loc_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     
        end
end
%% MIXED

%%
j=0;

eye_trace_nondef=nan([size(newdir_sel,1) 12]); 
eye_trace_alltrials=nan([size(newdir_sel,1) 12]);
for iAn=1:size(newdir_sel,2)
        j=j+1
        clear correct_value
        clear order
        frameRate=round(expt2_sel{iAn}.frameRate);
        epochs_length=length(resp_nondef_eye_trials{iAn});
        stims_value=stimulus{iAn}(1:end-1);
        correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
        for i=1:12;
        order(i)= find(stims_value==correct_value(i));
        end

        if epochs_length>0

eye_trace_nondef(iAn,:)=squeeze(nanmean(resp_nondef_eye_trials{iAn}(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(resp_nondef_eye_trials{iAn}(frameRate:2*frameRate,order)));
clear tmp_data
tmp_data=squeeze(nanmean(trials_eye_stims{iAn},3));
eye_trace_alltrials(iAn,:)=squeeze(nanmean(tmp_data(epochs_length-frameRate:end,order)))-...
squeeze(nanmean(tmp_data(frameRate:2*frameRate,order)));     

        end
end

%% LOPCOMOTION VS STILL TRIALS
figure,
clear tt
clear list_expt
for iAn=1:size(newdir_sel,2); tt(iAn)= sum(eye_trace_still(iAn,:)); end
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
for iAn=1:size(newdir_sel,2); tt(iAn)= sum(eye_trace_locomotion(iAn,:)); end
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
for iAn=1:size(newdir_sel,2); tt(iAn)= sum(eye_trace_nondef(iAn,:)); end
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
all_array_trials=nan(size(newdir_sel,2),15,12);
size(all_array_trials)

for iAn=1:size(newdir_sel,2)
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
   print(gcf,'-dpdf',[filepathanalysis, 'mergedplot_deltapupil_trials_average_40expt.pdf']);
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

   %% plot heatmaps all sessions
   data_trials=squeeze(nanmean(all_array_trials(:,:,:),1));
   size(data_trials)

    figure
    subplot(221)
    imagesc(data_trials(1:10,:),[-0.1 0.1]);
    colormap(bluered)
    colorbar
    hold on
    axis tight
    axis square 

%% plot same animals SEM
addpath 'G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals'
pattern_string='KS'
[animals_array]=find_sameanimals(newdir,pattern_string)
[x y]=find(animals_array==iAn)
%animals_array(iAn,~isnan(animals_array(iAn,:)))

for iAn=1:40
    [x y]=find(animals_array==iAn);
    x_animal(iAn)=x
    y_session(iAn)=y
end

clear data_trials_sameanimal
for n_animals=1:size(animals_array,1)
data_trials_sameanimal(n_animals,:,:)=...
    nanmean(all_array_trials(find(x_animal(:)==n_animals)',:,:),1);
end

size(data_trials_sameanimal)

%% plot for the same animals
   figure
    data_trials_sameanimals=squeeze(nanmean(data_trials_sameanimal(:,:,:),1));

    subplot(222)
    imagesc(data_trials_sameanimals(1:10,:),[-0.2 0.2]);
    colormap(bluered)
    colorbar
    hold on
    axis tight
    axis square 
    ylabel('#Trial');
    xlabel('Stimulus (deg)');
    title(sprintf('Delta pupil (mm^2) \n 13 mice'))
    set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(gca,'ticklength')*4);

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_diameter_trials\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'heatmap_deltapupil_av_trials_13mice_different_scale.pdf']);

%% statistical test
data_trials_sameanimal

fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');

axes(ax_heatmap_allmice_plot)
imagesc(squeeze(nanmean(data_trials_sameanimal(:,1:10,:),1)))
axis tight
axis square
cmap=colormap(bluered)
set(gca,'CLim',[-0.4 0.4])

cb=colorbar;
cb.Position = [.27,.082,.015,.05] 
caxis([-0.2 0.2]);

clear average_trials_animal

average_trials_animal=squeeze(nanmean(data_trials_sameanimal(:,1:10,:),1));


order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

clear av_nasal_trials
clear av_temporal_trials
   av_nasal_trials(:,:)= average_trials_animal(:,order_nasal);
   av_temporal_trials(:,:)= average_trials_animal(:,order_temporal);   

clear h0
clear pval
for ntrials=1:10;
    
    clear tmp_temporal
    clear tmp_nasal
    
tmp_temporal=av_temporal_trials(ntrials,:);
tmp_nasal=av_nasal_trials(ntrials,:);
[h0(ntrials), pval(ntrials)]=ttest2(tmp_nasal, tmp_temporal);
%[h0(ntrials), pval(ntrials)]=kstest2(tmp_nasal, tmp_temporal);

end

[h0; pval]
for ntrials=1:10
if pval(ntrials)<0.05 & h0(ntrials)==1
    str_name{ntrials}=sprintf('* %d', ntrials);    
else
    str_name{ntrials}=sprintf('  %d', ntrials);    
end
end
xtickstr={'0',[],[],'90',[],[],'180',[],[],'270',[],[]};
set(ax_heatmap_allmice_plot,'xtick',[1:1:12],'xticklabel',xtickstr,'ytick',[1:1:10],'yticklabel',str_name,...
     'tickdir','out','box','off','layer','top','color','none',...
     'fontsize',10,'ticklength',get(ax_heatmap_allmice_plot,'ticklength')*4);

set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals\plot_heatmaps_pupil_median_pvals\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_deltapupil_pupil_pvals_trials.pdf']);

