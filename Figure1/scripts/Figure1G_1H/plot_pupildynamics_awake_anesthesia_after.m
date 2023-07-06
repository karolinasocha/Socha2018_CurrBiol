%% PREPARE DATA
addpath(genpath('G:\mousebox\code\mouselab\users\karolina'))
eye_results_path='H:\local\users\karolina\170110_KS173_2P_KS_eye';
newdir_awake='170110_KS173_2P_KS\run03_ori12_V1_awake';
newdir_awake_after='170110_KS173_2P_KS\run03_ori12_V1_awake_after';
newdir_anesthesia='170110_KS173_2P_KS\run03_ori12_V1_anesthesia';

%%
expt_awake = frGetExpt(newdir_awake);
expt2_awake =doLoadStimLogs3(expt_awake);

expt_awake_after = frGetExpt(newdir_awake_after);
expt2_awake_after =doLoadStimLogs3(expt_awake_after);

expt_anesthesia = frGetExpt(newdir_anesthesia);
expt2_anesthesia =doLoadStimLogs3(expt_anesthesia);

%%
clear diameter_data_revision
filepathanalysis=['G:\mousebox\analysis\' newdir_anesthesia]
pupil_finalName = [filepathanalysis, '\diameter_data_revision.mat'];
load(pupil_finalName)
pupildiam_anesthesia=zscore(diameter_data_revision.diameter_no_artifact_smooth);

stimlog = expt2_anesthesia;
dtime = linspace(0,1,length(pupildiam_anesthesia));
p2time =  linspace(0,1,stimlog.nFramesTotal);
diam_final_anesthesia = pi*(interp1(dtime,pupildiam_anesthesia,p2time)/2).^2; % interpolated and changed mm square
%%
clear diameter_data_revision
filepathanalysis=['G:\mousebox\analysis\' newdir_awake]
pupil_finalName = [filepathanalysis, '\diameter_data_revision.mat'];
load(pupil_finalName)
pupildiam_awake=diameter_data_revision.diameter_no_artifact_smooth;
clear stimlog
stimlog = expt2_awake;
dtime = linspace(0,1,length(pupildiam_awake));
p2time =  linspace(0,1,stimlog.nFramesTotal);
diam_final_awake = pi*(interp1(dtime,pupildiam_awake,p2time)/2).^2; % interpolated and changed mm square

%%
clear diameter_data_revision
filepathanalysis=['G:\mousebox\analysis\' newdir_awake_after]
pupil_finalName = [filepathanalysis, '\diameter_data_revision.mat'];
load(pupil_finalName)
pupildiam_awake_after=diameter_data_revision.diameter_no_artifact_smooth;
clear stimlog
stimlog = expt2_awake_after;
dtime = linspace(0,1,length(pupildiam_awake_after));
p2time =  linspace(0,1,stimlog.nFramesTotal);
diam_final_awake_after = pi*(interp1(dtime,pupildiam_awake_after,p2time)/2).^2; % interpolated and changed mm square

%% epochs
[ddat_awake_after ddat_awake_after_stim]=calc_pupilepochs(expt2_awake_after,diam_final_awake_after);
[ddat_anesthesia ddat_anesthesia_stim] =calc_pupilepochs(expt2_anesthesia,diam_final_anesthesia);
[ddat_awake ddat_awake_stim]=calc_pupilepochs(expt2_awake,diam_final_awake);

%ddat_anesthesia_stimulus=tcEpochAverage2(diam_final_anesthesia', expt2_anesthesia.frames.stims)
%ddat_postawake_stimulus=tcEpochAverage2(diam_final_awake_after', expt2_awake_after.frames.stims)

%%
stimulus_values=[expt2_awake.prot.pars.ori];
stimulus_values=stimulus_values(1:end-1);
iStim_nasal=find(stimulus_values==0);
iStim_temporary=find(stimulus_values==180);
iStim_up=find(stimulus_values==90);
iStim_down=find(stimulus_values==270);
pre_frames=round(2*expt2_awake.frameRate)
post_frames=round(8*expt2_awake.frameRate)
stim_frames=pre_frames+round(5*expt2_awake.frameRate)
nasal_color=[0 0.5 1];
temporal_color=[0.5 0.5 0.5];
%%
fig(2) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
clf
ax_pupildynamics_allsessions_plot_nasal_awake = axes('position',[.15,.15,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_nasal_awake_after = axes('position',[.45,.15,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_nasal_anesthesia = axes('position',[.75,.15,.15,.27],'units','normalized');

%
ylimitation=[0 1.4];
tmp1=ddat_awake_after;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_nasal);
tmp_toplot_contrary=tmp1(:,:,iStim_temporary);

axes(ax_pupildynamics_allsessions_plot_nasal_awake_after)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight
set(ax_pupildynamics_allsessions_plot_nasal_awake_after,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_nasal_awake_after,'ticklength')*4);
ylim(ylimitation);
ylabel ('Awake After (mm^2)')
xlabel ('Time(s)')
title '0-180'
xlim([0 pre_frames+post_frames])

%
tmp1=ddat_awake;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_nasal);
tmp_toplot_contrary=tmp1(:,:,iStim_temporary);

axes(ax_pupildynamics_allsessions_plot_nasal_awake)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight
set(ax_pupildynamics_allsessions_plot_nasal_awake,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_nasal_awake,'ticklength')*4);
ylim(ylimitation);
ylabel ('Awake (mm^2)')
xlabel ('Time(s)')
title '0-180'
xlim([0 pre_frames+post_frames])
%
tmp1=ddat_anesthesia;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_nasal);
tmp_toplot_contrary=tmp1(:,:,iStim_temporary);

axes(ax_pupildynamics_allsessions_plot_nasal_anesthesia)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight;
hold on;
plot([pre_frames pre_frames],ylimitation,'--k');
hold on;
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight;
set(ax_pupildynamics_allsessions_plot_nasal_anesthesia,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_nasal_anesthesia,'ticklength')*4);
ylim(ylimitation);
ylabel ('Anesthesia (mm^2)');
xlabel ('Time(s)');
title '0-180';
xlim([0 pre_frames+post_frames]);
%% UPWARD DOWNWARD
% fig(2) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
%     'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);

ax_pupildynamics_allsessions_plot_up_awake = axes('position',[.15,.62,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_up_awake_after = axes('position',[.45,.62,.15,.27],'units','normalized');
ax_pupildynamics_allsessions_plot_up_anesthesia = axes('position',[.75,.62,.15,.27],'units','normalized');

%
ylimitation=[0 1.4];
tmp1=ddat_awake_after;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_up);
tmp_toplot_contrary=tmp1(:,:,iStim_down);

axes(ax_pupildynamics_allsessions_plot_up_awake_after)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight
set(ax_pupildynamics_allsessions_plot_up_awake_after,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_up_awake_after,'ticklength')*4);
ylim(ylimitation);
ylabel ('Awake After (mm^2)')
xlabel ('Time(s)')
title '90-270'
xlim([0 pre_frames+post_frames])

%
tmp1=ddat_awake;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_up);
tmp_toplot_contrary=tmp1(:,:,iStim_down);

axes(ax_pupildynamics_allsessions_plot_up_awake)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight
set(ax_pupildynamics_allsessions_plot_up_awake,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_up_awake,'ticklength')*4);
ylim(ylimitation);
ylabel ('Awake (mm^2)')
xlabel ('Time(s)')
title '90-270'
xlim([0 pre_frames+post_frames])
%
tmp1=ddat_anesthesia;
size(tmp1)
%
% NASAL-TEMPORAL -------------------------------------------------------------------------
%tmp1=squeeze(array_data(:,:));
tmp_toplot=tmp1(:,:,iStim_up);
tmp_toplot_contrary=tmp1(:,:,iStim_down);

axes(ax_pupildynamics_allsessions_plot_up_anesthesia)
hold on
h1=shadedErrorBar([1:length(tmp_toplot(:,:))],...
    nanmean(tmp_toplot(:,:),2),...
    nanstd(tmp_toplot(:,:),[],2)./...
    sqrt(size(tmp1,2)),'k');
hold on
h1.patch.FaceColor=nasal_color;
h1.patch.FaceAlpha=1;
hold on
h=shadedErrorBar([1:length(tmp_toplot_contrary(:,:))],...
    nanmean(tmp_toplot_contrary(:,:),2),...
    nanstd(tmp_toplot_contrary(:,:),[],2)./sqrt(size(tmp1,2)),'k')
h.patch.FaceColor=temporal_color;
h.patch.FaceAlpha=1;

axis tight
hold on
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
val=[-2,0,5,8];
axis tight
set(ax_pupildynamics_allsessions_plot_up_anesthesia,'ytick',[min(ylimitation):0.2:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_up_anesthesia,'ticklength')*4);
ylim(ylimitation);
ylabel ('Anesthesia (mm^2)')
xlabel ('Time(s)')
title '90-270'
xlim([0 pre_frames+post_frames])

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\Figures\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'example_pupil_awake_anesthesia_after.pdf']);


