%% Figure SM3 
% anesthesia recordings and post recordings

% clear all
% newdir_diameter_anesthesia={'181108_KS600_ephys_behav_JC\run01_direction_12dir_4Hz_ansth',...
%     '181108_KS602_ephys_behav_JC\run01_direction_12dir_4Hz_anesth'};
% 
% newdir_diameter_postawake={'181108_KS600_ephys_behav_JC\run02_direction_12dir_4Hz',...
%     '181108_KS602_ephys_behav_JC\run02_direction_12dir_4Hz'};

%%

clear all

newdir_diameter_anesthesia={'181108_KS600_ephys_behav_JC\run01_direction_12dir_4Hz_ansth',...
    '181108_KS602_ephys_behav_JC\run01_direction_12dir_4Hz_anesth',...
    '170110_KS173_2P_KS\run03_ori12_V1_anesthesia'};

newdir_diameter_postawake={'181108_KS600_ephys_behav_JC\run02_direction_12dir_4Hz',...
     '181108_KS602_ephys_behav_JC\run02_direction_12dir_4Hz',...
     '170110_KS173_2P_KS\run03_ori12_V1_awake_after'};
 
analysis_path='J:\data\analysis\';

for iAn=1:size(newdir_diameter_anesthesia,2)
   clear tmp_behav
   
tmp_behav=load([analysis_path newdir_diameter_anesthesia{iAn} '\behav_data.mat']);
behav{iAn}=tmp_behav.behav_data;
diameter_anesthesia_stim{iAn}=pi*(behav{iAn}.diameter.diam_stim./2).^2;
diameter_anesthesia_epochs{iAn}=pi*(behav{iAn}.diameter.diam_epochs./2).^2;
diameter_anesthesia{iAn}=behav{iAn}.diameter; % this is previous calculation in mm

end

for iAn=1:size(newdir_diameter_postawake,2)
   clear tmp_behav 
   
tmp_behav=load([analysis_path newdir_diameter_postawake{iAn} '\behav_data.mat']);
behav_awake{iAn}=tmp_behav.behav_data;
diameter_awake_stim{iAn}=pi*(behav_awake{iAn}.diameter.diam_stim./2).^2;
diameter_awake_epochs{iAn}=pi*(behav_awake{iAn}.diameter.diam_epochs./2).^2;
diameter_awake{iAn}=behav_awake{iAn}.diameter; % this is in mm

end
%% plot example

fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[30,30],'paperposition',[0,0,30,30]);
%
clf
ax_diameter_example = axes('position',[.2,.7,.4,.1],'units','normalized');

diam_ansth_example=pi*(diameter_anesthesia{1}.diameter./2).^2;
diam_awake_example=pi*(diameter_awake{1}.diameter./2).^2;

axes(ax_diameter_example)

plot(diam_awake_example./2,'b');
hold on
plot(diam_ansth_example./2,'k');
axis tight
ylim([0 0.8]);
hold on

plot([10+39265+round(100*39265/1.3234e+03) 10+39265],[0 0],'k','linewidth',3);
hold on
plot([10+39265 10+39265],[0 0.25],'k','linewidth',3);
title('Anesthesia Awake Pupil');
set(ax_diameter_example,'box','off','tickdir','out','ytick',[0:0.2:0.8],...
    'yticklabel',[0:0.2:0.8],...
    'ticklength',get(ax_diameter_example,'ticklength')*3); 

% to define xlims

%% PLO DELTA PUPIL FOR ANESTHESIA

eye_trace_awake_stims=nan([size(newdir_diameter_anesthesia,2)+1 12]);
eye_trace_anesthesia_stims=nan([size(newdir_diameter_anesthesia,2)+1 12]);

for iAn=1:size(newdir_diameter_anesthesia,2)
        clear temp
        frameRate=diameter_anesthesia{1}.cframespre/2;
        epochs_length=length(diameter_anesthesia_stim{iAn});

        temp=squeeze(nanmean(diameter_anesthesia_stim{iAn}(:,:,1:end-1),2));

        eye_trace_anesthesia_stims(iAn,:)=squeeze(nanmean(temp(epochs_length-frameRate:end,:)))-...
        squeeze(nanmean(temp(frameRate:2*frameRate,:)));
   
  clear temp
  try
        temp=squeeze(nanmean(diameter_awake_stim{iAn}(:,:,1:end-1),2));

        eye_trace_awake_stims(iAn,:)=squeeze(nanmean(temp(epochs_length-frameRate:end,:)))-...
        squeeze(nanmean(temp(frameRate:2*frameRate,:)));
  end
end
%%

%% get from the script G:\mousebox\code\mouselab\users\karolina\Socha2018_revision
% !!!! plot_pupildynamics_awake_anesthesia_after
% KS173 was with different format of data
ddat_anesthesia_stimulus=squeeze(nanmean(ddat_anesthesia_stim,2));
ddat_postawake_stimulus=squeeze(nanmean(ddat_awake_after_stim,2));
epochs_length=size(ddat_anesthesia_stimulus,1)
frameRate=round(expt2_anesthesia.frameRate);
stimulus_after=[expt2_anesthesia.prot.pars.ori];
[val ord]=sort(stimulus_after(1:end-1));

eye_trace_anesthesia_stims(size(eye_trace_anesthesia_stims,1),:)=squeeze(nanmean(ddat_anesthesia_stimulus(epochs_length-frameRate:end,ord)))-...
        squeeze(nanmean(ddat_anesthesia_stimulus(frameRate:2*frameRate,ord)));

epochs_length=size(ddat_postawake_stimulus,1)
frameRate=round(expt2_awake_after.frameRate);
stimulus_after=[expt2_awake_after.prot.pars.ori];
[val ord]=sort(stimulus_after(1:end-1));

eye_trace_awake_stims(size(eye_trace_awake_stims,1),:)=squeeze(nanmean(ddat_postawake_stimulus(epochs_length-frameRate:end,ord)))-...
        squeeze(nanmean(ddat_postawake_stimulus(frameRate:2*frameRate,ord)));

%%
figure
ax_diameter_delta = axes('position',[.2,.3,.2,.15],'units','normalized');

errorbar(1:size(eye_trace_awake_stims,2),squeeze(nanmean(eye_trace_awake_stims,1)),...
    nanstd(eye_trace_awake_stims,[],1)./sqrt(size(eye_trace_awake_stims,1)),'b');
hold on
plot([1 12],[0 0],'--k')
hold on
errorbar(1:size(eye_trace_anesthesia_stims,2),squeeze(nanmean(eye_trace_anesthesia_stims,1)),...
    nanstd(eye_trace_anesthesia_stims,[],1)./sqrt(size(eye_trace_anesthesia_stims,1)),'k');

axis tight

ylabel('Diff Pupil area all trials(mm^2)');
set(ax_diameter_delta,'ytick',[0:0.05:0.25],'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax_diameter_delta,'ticklength')*4);
ylim([-0.02 0.2]);
hold on

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'Figure1G_deltapupil_3mice_before_after_anesthesia_correct_final_figure_correct.pdf']);

%% pupil dynamics during anesthesia

%epochs
stimulus_after=[expt2_anesthesia.prot.pars.ori];
[val ord]=sort(stimulus_after(1:end-1));
trials_X_epochs_ordered{3}=ddat_anesthesia(:,:,[ord 13]);

trials_X_epochs_ordered{1}=diameter_anesthesia_epochs{1}(:,:,:);
trials_X_epochs_ordered{2}=diameter_anesthesia_epochs{2}(:,:,:);

average_diam_epochs_ordered{1}=squeeze(nanmean(trials_X_epochs_ordered{1},2));
average_diam_epochs_ordered{2}=squeeze(nanmean(trials_X_epochs_ordered{2},2));
average_diam_epochs_ordered{3}=squeeze(nanmean(trials_X_epochs_ordered{3},2));


for iAn=1:size(average_diam_epochs_ordered,2)
clear tmp1
tmp1=average_diam_epochs_ordered{iAn};
tmp=nan(1,4147);
    clear tempo
    tempo=squeeze(tmp1(:));
    dtime = linspace(0,1,length(tempo));
    p2time =  linspace(0,1,4147);
    tmp= interp1(dtime,tempo,p2time,'linear','extrap');
%tmp(:,iStim)=interp1(1:length(squeeze(tmp1(:,iStim))),squeeze(tmp1(:,iStim)),p2time);
%end
array_X_data(iAn,:)=tmp;
end

%%

expt2{1}=expt2_anesthesia
example_animal=1;
ylimitation=[0 1];
iStim_nasal=expt2{example_animal}.frames.stims{1,1}(1); % 0 deg
iStim_temporary=expt2{example_animal}.frames.stims{1,7}(1); % 180 deg
iStim_up=expt2{example_animal}.frames.stims{1,4}(1); % 90 deg
iStim_down=expt2{example_animal}.frames.stims{1,10}(1); % 270 deg
iStim_300=expt2{example_animal}.frames.stims{1,11}(1); % 300 deg
iStim_330=expt2{example_animal}.frames.stims{1,12}(1); % 330 deg
iStim_30=expt2{example_animal}.frames.stims{1,2}(1); % 30 deg
iStim_60=expt2{example_animal}.frames.stims{1,3}(1); % 60 deg
iStim_120=expt2{example_animal}.frames.stims{1,5}(1); % 120 deg
iStim_150=expt2{example_animal}.frames.stims{1,6}(1); % 150 deg
iStim_210=expt2{example_animal}.frames.stims{1,8}(1); % 210 deg
iStim_240=expt2{example_animal}.frames.stims{1,9}(1); % 240 deg

pre_frames=round(2*expt2{example_animal}.frameRate)
post_frames=round(8*expt2{example_animal}.frameRate)
stim_frames=pre_frames+round(5*expt2{example_animal}.frameRate)
nasal_color=[0 0.5 1];
temporal_color=[0.5 0.5 0.5];
%
fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
ax_pupildynamics_allsessions_plot_nasal = axes('position',[.1,.07,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_up = axes('position',[.30,.07,.15,.27],'units','normalized');

ax_pupildynamics_allsessions_plot_300 = axes('position',[.1,.37,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_330 = axes('position',[.30,.37,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_30 = axes('position',[.1,.67,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_60 = axes('position',[.30,.67,.15,.27],'units','normalized');

% axes(ax_pupildynamics_allsessions_plot_nasal)
% pupildynamic_mean_nasal=squeeze(nanmean(array_data(:,:),1));
%tmp1=squeeze(array_data(:,:));
%tmp1=array_X_data_sameanimals;
array_X_data_sameanimals=array_X_data;
tmp1=array_X_data
nanimals=size(array_X_data,1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,iStim_nasal-pre_frames:iStim_nasal+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_temporary-pre_frames:iStim_temporary+post_frames);

axes(ax_pupildynamics_allsessions_plot_nasal)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_nasal,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_nasal,'ticklength')*4);
ylim(ylimitation);
ylabel ('Pupil (norm)')
xlabel ('Time(s)')
title '0-180'
xlim([0 pre_frames+post_frames])
% UP-DOWN -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_up-pre_frames:iStim_up+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_down-pre_frames:iStim_down+post_frames);

axes(ax_pupildynamics_allsessions_plot_up)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_up,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_up,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '90-270';
%
% other directions
% 300-120 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_300-pre_frames:iStim_300+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_120-pre_frames:iStim_120+post_frames);

axes(ax_pupildynamics_allsessions_plot_300)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_300,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_300,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '300-120';
% 
% 330-150 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_330-pre_frames:iStim_330+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_150-pre_frames:iStim_150+post_frames);

axes(ax_pupildynamics_allsessions_plot_330)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_330,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_330,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '330-150';

% 30-210 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_30-pre_frames:iStim_30+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_210-pre_frames:iStim_210+post_frames);

axes(ax_pupildynamics_allsessions_plot_30)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_30,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_30,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '30-210';

% 60-240 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_60-pre_frames:iStim_60+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_240-pre_frames:iStim_240+post_frames);

axes(ax_pupildynamics_allsessions_plot_60)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_60,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_60,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '60-240';

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_anesthesia_postAwake\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'summary_pupildynamic_3mice_anesthesia_mm2.pdf']);

%% post awake

%% pupil dynamics during anesthesia

%epochs
stimulus_after=[expt2_awake_after.prot.pars.ori];
[val ord]=sort(stimulus_after(1:end-1));
trials_X_epochs_ordered_awake{2}=ddat_awake_after(:,:,[ord 13]);

trials_X_epochs_ordered_awake{1}=diameter_awake_epochs{1}(:,:,:);
%trials_X_epochs_ordered_awake{2}=diameter_anesthesia_epochs{2}(:,:,:);

average_diam_epochs_ordered_awake{1}=squeeze(nanmean(trials_X_epochs_ordered_awake{1},2));
average_diam_epochs_ordered_awake{2}=squeeze(nanmean(trials_X_epochs_ordered_awake{2},2));
%average_diam_epochs_ordered{3}=squeeze(nanmean(trials_X_epochs_ordered_awake{3},2));


for iAn=1:size(average_diam_epochs_ordered_awake,2)
clear tmp1
tmp1=average_diam_epochs_ordered_awake{iAn};
tmp=nan(1,4147);
    clear tempo
    tempo=squeeze(tmp1(:));
    dtime = linspace(0,1,length(tempo));
    p2time =  linspace(0,1,4147);
    tmp= interp1(dtime,tempo,p2time,'linear','extrap');
%tmp(:,iStim)=interp1(1:length(squeeze(tmp1(:,iStim))),squeeze(tmp1(:,iStim)),p2time);
%end
array_X_data_awake(iAn,:)=tmp;
end

%%
example_animal=1;
ylimitation=[0 1.5];
iStim_nasal=expt2{example_animal}.frames.stims{1,1}(1); % 0 deg
iStim_temporary=expt2{example_animal}.frames.stims{1,7}(1); % 180 deg
iStim_up=expt2{example_animal}.frames.stims{1,4}(1); % 90 deg
iStim_down=expt2{example_animal}.frames.stims{1,10}(1); % 270 deg
iStim_300=expt2{example_animal}.frames.stims{1,11}(1); % 300 deg
iStim_330=expt2{example_animal}.frames.stims{1,12}(1); % 330 deg
iStim_30=expt2{example_animal}.frames.stims{1,2}(1); % 30 deg
iStim_60=expt2{example_animal}.frames.stims{1,3}(1); % 60 deg
iStim_120=expt2{example_animal}.frames.stims{1,5}(1); % 120 deg
iStim_150=expt2{example_animal}.frames.stims{1,6}(1); % 150 deg
iStim_210=expt2{example_animal}.frames.stims{1,8}(1); % 210 deg
iStim_240=expt2{example_animal}.frames.stims{1,9}(1); % 240 deg

pre_frames=round(2*expt2{example_animal}.frameRate)
post_frames=round(8*expt2{example_animal}.frameRate)
stim_frames=pre_frames+round(5*expt2{example_animal}.frameRate)
nasal_color=[0 0.5 1];
temporal_color=[0.5 0.5 0.5];
%
fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
ax_pupildynamics_allsessions_plot_nasal = axes('position',[.1,.07,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_up = axes('position',[.30,.07,.15,.27],'units','normalized');

ax_pupildynamics_allsessions_plot_300 = axes('position',[.1,.37,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_330 = axes('position',[.30,.37,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_30 = axes('position',[.1,.67,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_60 = axes('position',[.30,.67,.15,.27],'units','normalized');

% axes(ax_pupildynamics_allsessions_plot_nasal)
% pupildynamic_mean_nasal=squeeze(nanmean(array_data(:,:),1));
%tmp1=squeeze(array_data(:,:));
%tmp1=array_X_data_sameanimals;
array_X_data_sameanimals=array_X_data_awake;
tmp1=array_X_data_awake
nanimals=size(array_X_data_awake,1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,iStim_nasal-pre_frames:iStim_nasal+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_temporary-pre_frames:iStim_temporary+post_frames);

axes(ax_pupildynamics_allsessions_plot_nasal)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_nasal,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_nasal,'ticklength')*4);
ylim(ylimitation);
ylabel ('Pupil (norm)')
xlabel ('Time(s)')
title '0-180'
xlim([0 pre_frames+post_frames])
% UP-DOWN -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_up-pre_frames:iStim_up+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_down-pre_frames:iStim_down+post_frames);

axes(ax_pupildynamics_allsessions_plot_up)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_up,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_up,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '90-270';
%
% other directions
% 300-120 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_300-pre_frames:iStim_300+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_120-pre_frames:iStim_120+post_frames);

axes(ax_pupildynamics_allsessions_plot_300)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_300,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_300,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '300-120';
% 
% 330-150 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_330-pre_frames:iStim_330+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_150-pre_frames:iStim_150+post_frames);

axes(ax_pupildynamics_allsessions_plot_330)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_330,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_330,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '330-150';

% 30-210 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_30-pre_frames:iStim_30+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_210-pre_frames:iStim_210+post_frames);

axes(ax_pupildynamics_allsessions_plot_30)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_30,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_30,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '30-210';

% 60-240 -------------------------------------------------------------------------
tmp_toplot=tmp1(:,iStim_60-pre_frames:iStim_60+post_frames);
tmp_toplot_contrary=tmp1(:,iStim_240-pre_frames:iStim_240+post_frames);

axes(ax_pupildynamics_allsessions_plot_60)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),1),...
    nanstd(tmp_toplot(:,:),[],1)./...
    sqrt(nanimals),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),1),nanstd(tmp_toplot_contrary(:,:),[],1)./...
    sqrt(nanimals),'k');
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
set(ax_pupildynamics_allsessions_plot_60,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_60,'ticklength')*4);
ylim(ylimitation);
%ylabel ('Pupil (norm)');
xlabel ('Time(s)');
title '60-240';

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\pupil_anesthesia_postAwake\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'summary_pupildynamic_3mice_anesthesia_mm2.pdf']);
