%% run script load_data_cortical_fractionrunning.m to load all dataset
% it will average population reponse over all experiments

%%
clear all
file_path='G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\dFF_cortical_data_fractionrunning_heatmaps'
newdir='180822_KS501_2P_JC\run00_direction_12dir_4Hz'

%%
clear timecourses
timecourses=load(['J:\data\analysis\' newdir '\timecourses_02.mat']);
timecourses=timecourses.timecourses_02;
for iPlane=1:4
tcs.epochs{iPlane}=timecourses.calcium{iPlane}.epoch;
tcs.stims{iPlane}=timecourses.calcium{iPlane}.stim;
end

%%

clear timecourses
timecourses=load(['J:\data\analysis\' newdir '\timecourses.mat']);
for iPlane=1:4
clear responsive_elements
responsive_elements=timecourses.timecourses.boutons_selection.responsive_elements{iPlane}
tcs.epochs{iPlane}=tcs.epochs{iPlane}(:,:,:,responsive_elements);
tcs.stims{iPlane}=tcs.stims{iPlane}(:,:,:,responsive_elements);
end

w(1)=0;
ws(1)=0;
clear tmp_tc
for iPlane=1:4
[x1,y1,z1,w(1+iPlane)]=size(tcs.epochs{iPlane});
[x1s,y1s,z1s,ws(1+iPlane)]=size(tcs.stims{iPlane});
end
tmp_tc=nan(x1,y1,z1,sum(w));
tmp_tc_stims=nan(x1s,y1s,z1s,sum(ws));
for iPlane=1:4
    iiplane=iPlane+1
tmp_tc(:,:,:,1+sum(w(1:iPlane)):sum(w(1:iiplane)))=tcs.epochs{iPlane};
tmp_tc_stims(:,:,:,1+sum(w(1:iPlane)):sum(w(1:iiplane)))=tcs.stims{iPlane};

tc_stim_allplanes_stims=tmp_tc_stims;
tc_stim_allplanes_epochs=tmp_tc;
end

average_response=squeeze(nanmean(tc_stim_allplanes_epochs,2));
population_response_prctile50=squeeze(prctile(average_response,50,3));
population_response_prctile75=squeeze(prctile(average_response,75,3));
population_response_prctile5=squeeze(prctile(average_response,5,3));
population_response_prctile25=squeeze(prctile(average_response,25,3));

%%
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
%
clf
for iStim=1:12
ax_population_responses_prctile75(iStim) = axes('position',[0.07*(iStim-1)+.1,.45,.05,.07],'units','normalized');
ax_population_responses_prctile5(iStim) = axes('position',[0.07*(iStim-1)+.1,.25,.05,.07],'units','normalized');
ax_population_responses_prctile50(iStim) = axes('position',[0.07*(iStim-1)+.1,.65,.05,.07],'units','normalized');
ax_population_responses_prctile25(iStim) = axes('position',[0.07*(iStim-1)+.1,.85,.05,.07],'units','normalized');

end

%
ylimitation_prctile5=[-5 10];
ylimitation_prctile75=[5 25];
ylimitation_prctile50=[0 25];
ylimitation_prctile25=[0 25];

for iStim=1:12
axes(ax_population_responses_prctile75(iStim))
plot(population_response_prctile75(:,iStim)','color','r','linewidth',1);
hold on
plot(population_response_prctile25(:,iStim)','color','b','linewidth',1);
hold on
plot(population_response_prctile50(:,iStim)','color','k','linewidth',1);
hold on
axis tight
set(ax_population_responses_prctile75(iStim),'ytick',[min(ylimitation_prctile75):5:max(ylimitation_prctile75)],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(ax_population_responses_prctile75(iStim),'ticklength')*4);
hold on
plot([16 16],ylimitation_prctile75,'--k');
hold on
plot([55 55],ylimitation_prctile75,'--k');
axis off
end

axes(ax_population_responses_prctile75(12));
plot([80 80],[5 10],'-k','linewidth',3);
text(80,10,'5% dF/F');

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\dFF_cortical_data_fractionrunning_heatmaps\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'example_fig6_KS5001_population_responses.pdf']);
