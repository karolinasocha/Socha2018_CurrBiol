% 
clear all
newdir_random={'180823_KS501_2P_JC\run01_direction_12dir_4Hz_random',...
 '180824_KS501_2P_JC\run01_direction_12dir_4Hz_random',...
 '180823_KS503_2P_JC\run01_direction_12dir_4Hz_random',...
 '180824_KS503_2P_JC\run01_direction_12dir_4Hz_random'};

analysis_path='J:\data\analysis\';

for iAn=1:size(newdir_random,2)
   clear tmp_behav 
tmp_behav=load([analysis_path newdir_random{iAn} '\behav_data.mat']);
behav{iAn}=tmp_behav.behav_data;
end

clear diam_animals
clear velo_animals
clear velo_epochs
clear diam_epochs
for iAn=1:size(newdir_random,2)
velo_animals{iAn}=behav{iAn}.velocity.vel_epochs;
diam_animals{iAn}=behav{iAn}.diameter.diam_epochs;
end

stimulus=mod([0:30:330],360);
%% interpolate

for iAn=1:size(newdir_random,2)
    clear dtime
    clear p2time
    dtime = linspace(0,1,size(diam_animals{iAn},1));
    p2time =  linspace(0,1,size(diam_animals{1},1));
for iStim = 1:13
    for iTrial=1:10
        clear tempo
        tempo=diam_animals{iAn}(:,iTrial,iStim);
        diam_interpolated{iAn}(:,iTrial,iStim)=interp1(dtime,tempo,p2time,'linear','extrap');
    end
end

end

diam_epochs=cell2mat(diam_interpolated);
diam_epochs_random=cell2mat(diam_interpolated);
velo_epochs=cell2mat(velo_animals);
%%
%% random experiments
for iAn=1:4
clear rand_diam
rand_diam=squeeze(nanmean(behav{iAn}.diameter.diam_stim,1));
rand_diam_nasal{iAn}=rand_diam(:,[1,2,3,4,11,12]);
rand_diam_tempo{iAn}=rand_diam(:,[5,6,7,8,9,10]);
end
% 
% save('G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\rand_diam_nasal.mat','rand_diam_nasal');
% save('G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\rand_diam_tempo.mat','rand_diam_tempo');

%% responses pupil
stimulus=[0,30,60,90,120,150,180,210,240,270,300,330];
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
%
clf

for iStim=1:6
ax_avpupil_stims(iStim) = axes('position',[.045+.165*(iStim-1),.075,.13,.13],'units','normalized');
ax_avpupil_stims(iStim+6) = axes('position',[.045+.165*(iStim-1),.3,.13,.13],'units','normalized');
end

%
ntrials=size(velo_epochs,2)
nstims=size(velo_epochs,3)

for iStim=1:nstims-1
axes(ax_avpupil_stims(iStim));
hold all
clear x
clear y
clear errBar
x=1:size(diam_epochs,1);
y=squeeze(diam_epochs(:,:,iStim));
errBar=nanstd(y,[],2)./sqrt(ntrials);
H=shadedErrorBar(x,nanmean(y,2),errBar,'r');
hold on
cframespre=behav{iAn}.diameter.cframespre;
cframesstim=cframespre+behav{iAn}.diameter.cframesstim;
plot([cframespre cframespre],[0.5 2.5],'--','color','k');
plot([cframesstim cframesstim],[0.5 2.5],'--','color','k');
axis tight;
ylim([0.8 1.2]);

title(sprintf('Stim %d',stimulus(iStim)))
set(ax_avpupil_stims(iStim),'box','off','tickdir','out','xtick',[cframespre cframesstim],...
   'xticklabel',[0 5],'ticklength',get(ax_avpupil_stims(iStim),'ticklength')*4);
end

%%
set(gcf,'paperunits','centimeters','papersize' ,[22,22],'color','w','paperposition',[0,0,21,21],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\']; 
print(gcf,'-dpdf',[filepathanalysis, '\pupil_stimulus_results_3mice_rand_4sessions.pdf']);

%% %
newdir_nonrandom={'180823_KS501_2P_JC\run00_direction_12dir_4Hz',...
 '180824_KS501_2P_JC\run00_direction_12dir_4Hz',...
 '180823_KS503_2P_JC\run00_direction_12dir_4Hz',...
 '180824_KS503_2P_JC\run00_direction_12dir_4Hz'};


analysis_path='J:\data\analysis\';

for iAn=1:size(newdir_nonrandom,2)
   clear tmp_behav 
tmp_behav=load([analysis_path newdir_nonrandom{iAn} '\behav_data.mat']);
behav{iAn}=tmp_behav.behav_data;
end

clear diam_animals
clear velo_animals
clear velo_epochs
clear diam_epochs
for iAn=1:size(newdir_nonrandom,2)
velo_animals{iAn}=behav{iAn}.velocity.vel_epochs;
diam_animals{iAn}=behav{iAn}.diameter.diam_epochs;
end

%% interpolate
clear diam_interpolated
for iAn=1:size(newdir_nonrandom,2)
    clear dtime
    clear p2time
    dtime = linspace(0,1,size(diam_animals{iAn},1));
    p2time =  linspace(0,1,size(diam_animals{3},1));
for iStim = 1:13
    for iTrial=1:10
        clear tempo
        tempo=diam_animals{iAn}(:,iTrial,iStim);
        size(tempo)
        diam_interpolated{iAn}(:,iTrial,iStim)=interp1(dtime,tempo,p2time,'linear','extrap');
    end
end

end

diam_epochs=cell2mat(diam_interpolated);
diam_epochs_non_random=cell2mat(diam_interpolated);
velo_epochs=cell2mat(velo_animals);
%%
%% nonrandom experiments
for iAn=1:4
clear nonrand_diam
nonrand_diam=squeeze(nanmean(behav{iAn}.diameter.diam_stim,1));
nonrand_diam_nasal{iAn}=nonrand_diam(:,[1,2,3,4,11,12]);
nonrand_diam_tempo{iAn}=nonrand_diam(:,[5,6,7,8,9,10]);
end

% save('G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\nonrand_diam_nasal.mat','nonrand_diam_nasal');
% save('G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\nonrand_diam_tempo.mat','nonrand_diam_tempo');

%% responses pupil
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
%
clf

for iStim=1:6
ax_avpupil_stims(iStim) = axes('position',[.045+.165*(iStim-1),.075,.13,.13],'units','normalized');
ax_avpupil_stims(iStim+6) = axes('position',[.045+.165*(iStim-1),.3,.13,.13],'units','normalized');
end

%
stimulus=mod([0:30:330],360);
ntrials=size(velo_epochs,2)
nstims=size(velo_epochs,3)

for iStim=1:nstims-1
axes(ax_avpupil_stims(iStim));
hold all
clear x
clear y
clear errBar
x=1:size(diam_epochs,1);
y=squeeze(diam_epochs(:,:,iStim));
errBar=nanstd(y,[],2)./sqrt(ntrials);
H=shadedErrorBar(x,nanmean(y,2),errBar,'k');
hold on
cframespre=behav{iAn}.diameter.cframespre;
cframesstim=cframespre+behav{iAn}.diameter.cframesstim;
plot([cframespre cframespre],[0.5 2.5],'--','color','k');
plot([cframesstim cframesstim],[0.5 2.5],'--','color','k');
axis tight;
ylim([0.8 1.2]);

title(sprintf('Stim %d',stimulus(iStim)))
set(ax_avpupil_stims(iStim),'box','off','tickdir','out','xtick',[cframespre cframesstim],...
   'xticklabel',[0 5],'ticklength',get(ax_avpupil_stims(iStim),'ticklength')*4);
end

%%
set(gcf,'paperunits','centimeters','papersize' ,[22,22],'color','w','paperposition',[0,0,21,21],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\2p_analysis\pipeline_randomstimulus_behavior\']; 
print(gcf,'-dpdf',[filepathanalysis, '\pupil_stimulus_results_3mice_nonrand_4sessions.pdf']);


%%

%% responses pupil
stimulus=[0,30,60,90,120,150,180,210,240,270,300,330];
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
%
clf

for iStim=1:6
ax_avpupil_stims(iStim) = axes('position',[.045+.165*(iStim-1),.075,.1,.15],'units','normalized');
ax_avpupil_stims(iStim+6) = axes('position',[.045+.165*(iStim-1),.3,.1,.15],'units','normalized');
end

%
ntrials=size(velo_epochs,2)
nstims=size(velo_epochs,3)

for iStim=1:nstims-1
axes(ax_avpupil_stims(iStim));
hold all
clear x
clear y
clear errBar
x=1:size(diam_epochs_non_random,1);
y=squeeze(diam_epochs_non_random(:,:,iStim));
errBar=nanstd(y,[],2)./sqrt(ntrials);
H=shadedErrorBar(x,nanmean(y,2),errBar,'k');
H.patch.FaceAlpha=0.3;
hold on

clear x
clear y
clear errBar
x=1:size(diam_epochs_random,1);
y=squeeze(diam_epochs_random(:,:,iStim));
errBar=nanstd(y,[],2)./sqrt(ntrials);
H2=shadedErrorBar(x,nanmean(y,2),errBar,'r');
H2.patch.FaceAlpha=0.3;
hold on
cframespre=behav{iAn}.diameter.cframespre;
cframesstim=cframespre+behav{iAn}.diameter.cframesstim;
plot([cframespre cframespre],[0.5 2.5],'--','color','k');
plot([cframesstim cframesstim],[0.5 2.5],'--','color','k');

axis tight;
ylim([0.7 1.2]);

title(sprintf('Stim %d',stimulus(iStim)))
set(ax_avpupil_stims(iStim),'box','off','tickdir','out','xtick',[cframespre cframesstim],...
   'xticklabel',[0 5],'ticklength',get(ax_avpupil_stims(iStim),'ticklength')*4);
end

set(gcf,'paperunits','centimeters','papersize' ,[22,22],'color','w','paperposition',[0,0,21,21],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\FigureS1C_random_nonrandom\']; 
print(gcf,'-dpdf',[filepathanalysis, '\pupil_stimulus_results_3mice_nonrand_random_4sessions.pdf']);
print(gcf,'-dpng',[filepathanalysis, '\pupil_stimulus_results_3mice_nonrand_random_4sessions.png']);
