%% zscoring anesthesia data to plot pupil dynamics in zscore  
%double check data set if correctly calculated - 190520
clear all
%newdir='190216_XH315_behav_JC\run01_direction_cardinal_4Hz_fullsize';

pathname_logfile='Z:\data\presentation\';
pathname_eye='Z:\data\eyecam\';
pathname_analysis='J:\data\analysis\'

newdir_diameter_anesthesia={'181108_KS600_ephys_behav_JC\run01_direction_12dir_4Hz_ansth',...
    '181108_KS602_ephys_behav_JC\run01_direction_12dir_4Hz_anesth'};

% newdir_diameter_postawake={'181108_KS600_ephys_behav_JC\run02_direction_12dir_4Hz',...
%     '181108_KS602_ephys_behav_JC\run02_direction_12dir_4Hz'};

analysis_path='J:\data\analysis\';


for iAn=1:size(newdir_diameter_anesthesia,2)
   clear tmp_behav
   
tmp_behav=load([analysis_path newdir_diameter_anesthesia{iAn} '\behav_data.mat']);
behav{iAn}=tmp_behav.behav_data;
X=behav{iAn}.diameter.diameter; % this is previous calculation in mm

Z{iAn} = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X)) , ...
                     nanstd(X));
end

clear log_file
for iAn=1:size(newdir_diameter_anesthesia,2)
   clear tmp_behav 
   clear log
 tmp_behav=load([analysis_path newdir_diameter_anesthesia{iAn} '\behav_data.mat']);
 behav{iAn}=tmp_behav.behav_data;
load([analysis_path newdir_diameter_anesthesia{iAn} '\log.mat']);
log_file{iAn}=log2
end

%%
for iAn=1:2;

import code2p.thirdparty.pyvstim.parseLog
log = log_file{iAn}
% compare with pupil_diameter_trace
camtime = interp1(log.cam1(:,4),log.cam1(:,3)/1000.,...
    1:length(Z{iAn}));

flips = log.screen(:,3)/1000.;
flipidx = log.screen(:,4);
stimtimes = [];
uflips = unique(flipidx);
for i = 1:length(uflips)
    stimtimes = [stimtimes,flips(flipidx == uflips(i))];
end
stimduration = max(ceil(diff(stimtimes,1)));
blanks = ~log.vstim(:,6);
oris = log.vstim(:,8);
iTrials = log.vstim(:,4);
iStims = log.vstim(:,3);
ptimes = log.vstim(:,2);
stim_size = log.vstim(:,12);

stimcode = iStims(find(diff(blanks)>0)+1) +1;
stimori = oris(find(diff(blanks)>0)+1);
trialcode = iTrials(find(diff(blanks)>0)+1);
stimdur = unique(round(diff(stimtimes,1)));

nstims = length(unique(stimcode));
trialidx = find(stimcode == 1);
ntrials = length(trialidx);

cframespre = behav{iAn}.diameter.cframespre; % 2 sec before
cframespost = behav{iAn}.diameter.cframespost; % 3 sec after ending stimulus
camidx = -1*cframespre:cframespost; % length of epoch
ddat = nan(length(camidx),ntrials,nstims);
size(ddat)

iinan = find(~isnan(Z{iAn}));
diameter=interp1(iinan,Z{iAn}(iinan),1:length(Z{iAn}));
%
for iStim = unique(stimcode)'
    trialidx = find(stimcode == iStim);
    ntrials = length(trialidx);
    for iTrial = 1:ntrials              
            
          icamFrame  = find(camtime <= stimtimes(1,trialidx(iTrial)),1,'last');
          if camidx(end)+icamFrame < length(diameter)
            ddat(:,iTrial,iStim) = diameter(icamFrame + camidx); 
          else
              disp(['Skipped trial ',num2str(iTrial),' ' ,num2str(iStim)])
          end
             hold all
           plot(ddat(:,iTrial,iStim))    
           
    end
end
%

pupildata.diameter.time=camtime;
pupildata.diameter.diameter=diameter;
pupildata.diameter.diam_epochs=ddat;
pupildata.diameter.diam_prestim=ddat(1:cframespre,:,:);
pupildata.diameter.info='2sec prestim; 5 sec stim; 3 sec post stim';

pupildata.diameter.cframespre=cframespre;
pupildata.diameter.cframespost=cframespost;
pupildata.stimulus=stimori;

pupildata.newdir=newdir_diameter_anesthesia{iAn};
pupildata_zscored=pupildata;
filepathanalysis=['J:\data\analysis\' newdir_diameter_anesthesia{iAn}]; 
pupil_name = [filepathanalysis, '\pupildata_zscored.mat'];
%save(pupil_name,'-v7.3','pupildata_zscored');
diameter_zscored{iAn}=ddat;
end

%% PREPARE DATA

eye_results_path='H:\local\users\karolina\170110_KS173_2P_KS_eye';
newdir_awake='170110_KS173_2P_KS\run03_ori12_V1_awake';
newdir_awake_after='170110_KS173_2P_KS\run03_ori12_V1_awake_after';
newdir_anesthesia='170110_KS173_2P_KS\run03_ori12_V1_anesthesia';

%
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
X=diameter_data_revision.diameter_no_artifact_smooth;
Z = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X)) , ...
                    nanstd(X));

iinan = find(~isnan(Z));
diameter_nan=interp1(iinan,Z(iinan),1:length(Z));

stimlog = expt2_anesthesia;
dtime = linspace(0,1,length(diameter_nan));
p2time =  linspace(0,1,stimlog.nFramesTotal);
diam_final_anesthesia = interp1(dtime,diameter_nan,p2time); % interpolated and changed mm square
[ddat_anesthesia ddat_anesthesia_stim] =calc_pupilepochs(expt2_anesthesia,diam_final_anesthesia);

ddat_anesthesia_stimulus=squeeze(nanmean(ddat_anesthesia,2));
epochs_length=size(ddat_anesthesia_stimulus,1)
frameRate=round(expt2_anesthesia.frameRate);
stimulus_after=[expt2_anesthesia.prot.pars.ori];
[val ord]=sort(stimulus_after(1:end-1));

ddat_anesthesia_ordered=ddat_anesthesia(:,:,[ord 13]);
%%

diameter_zscored{3}=ddat_anesthesia_ordered;

average_diam_epochs_ordered{1}=squeeze(nanmean(diameter_zscored{1},2));
average_diam_epochs_ordered{2}=squeeze(nanmean(diameter_zscored{2},2));
average_diam_epochs_ordered{3}=squeeze(nanmean(diameter_zscored{3},2));

for iAn=1:size(average_diam_epochs_ordered,2)
clear tmp1
clear tmp
tmp1=average_diam_epochs_ordered{iAn};
tmp=nan(315,13); 
for nStim=1:13
    clear tempo
    tempo=squeeze(tmp1(:,nStim));
    dtime = linspace(0,1,length(tempo));
    p2time =  linspace(0,1,315);
    tmp(:,nStim)= interp1(dtime,tempo,p2time,'linear','extrap');
end
array_X_data(iAn,:,:)=tmp;
end
%%

example=squeeze(nanmean(array_X_data,1));
figure;
plot(example(:,1))
axis tight
%% PLOT ZSCORED VALUE

fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
ax_pupildynamics_allsessions_plot_anesthesia = axes('position',[.15,.2,.06,.1],'units','normalized');
ylimitation=[-1 1.5];
%ylimitation=[0.4 0.8];

frameRate=315/10;

pre_frames=round(2*frameRate)
post_frames=round(8*frameRate)
stim_frames=pre_frames+round(5*frameRate)
val=[-2.2,0,5,8];

%tmp_toplot=array_X_data_sameanimals_stationary_nasal(:,iStim_nasal-pre_frames:iStim_nasal+post_frames);
axes(ax_pupildynamics_allsessions_plot_anesthesia)
clear tmp_nasal_av_anesthesia
clear tmp_nasal_std_anesthesia
tmp_nasal_av_anesthesia=nanmean(array_X_data(:,:,1),1);
tmp_nasal_std_anesthesia=nanstd(array_X_data(:,:,1),[],1)./...
    sqrt(3);

clear tmp_temporal_av_loc
clear tmp_temporal_std_loc

clear tmp_temporal_av_anesthesia
clear tmp_temporal_std_anesthesia
tmp_temporal_av_anesthesia=nanmean(array_X_data(:,:,7),1);
tmp_temporal_std_anesthesia=nanstd(array_X_data(:,:,7),[],1)./...
    sqrt(3);

hold on
h1=shadedErrorBar([1:length(tmp_nasal_av_anesthesia(:))],...
    tmp_nasal_av_anesthesia,...
    tmp_nasal_std_anesthesia,'b');

hold on

h2=shadedErrorBar([1:length(tmp_temporal_av_anesthesia(:))],...
    tmp_temporal_av_anesthesia,...
    tmp_temporal_std_anesthesia,'k');

hold on
axis tight
plot([pre_frames pre_frames],ylimitation,'--k');
hold on
plot([stim_frames stim_frames],ylimitation,'--k');
set(ax_pupildynamics_allsessions_plot_anesthesia,'ytick',[min(ylimitation):0.5:max(ylimitation)],'xtick',[0,pre_frames,stim_frames,pre_frames+post_frames],'xticklabel',val,...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_pupildynamics_allsessions_plot_anesthesia,'ticklength')*4);
title 'Anesthesia'
ylabel('Pupil (zscored)')
xlabel('Time')
%

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\plot_pupildynmics_locomotion_stationary\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'plot_pupildynamics_anesthesia_temporal_zscored.pdf']);

