% newdirs={'240314_AE01__no_cam_SH', '240305_MG04__no_cam_SH','240311_SH056__no_cam_SH'}
% sessions={'run00_140423_12dir-sf08-tf4SinSeq5s-20T, 'run00_195835_12dir-sf08-tf-4SinSeq''}
% Experiments list

% EFFECT IN FEMALE MICE
% '240314_AE01__no_cam_SH\run00_140423_12dir-sf08-tf4SinSeq5s-20T' % % -FEMALE; sequence; 10 trials; poor quality to remove
% '240305_MG04__no_cam_SH\run00_195835_12dir-sf08-tf-4SinSeq' % -FEMALE; % sequence; 10 trials; no effect

% EFFECT DIFFERENT DURATION
% '240311_SH056__no_cam_SH\run00_154140_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH056__no_cam_SH\run01_160616_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH056__no_cam_SH\run02_162601_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH056__no_cam_SH\run03_164335_12dir-sf08-tf4SinSeq2s' - 2s duration

% '240311_SH057__no_cam_SH\run00_170314_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH057__no_cam_SH\run01_172413_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH057__no_cam_SH\run02_174356_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH057__no_cam_SH\run03_180209_12dir-sf08-tf4SinSeq2s' - 3s duration

% EFFECT DIFFERENT TFs
% '240313_SH055_no_cam_SH\run00_133318_5tf-4dir-sf08-SinRand-5s' 
% '240313_SH056__no_cam_SH\run00_145132_5tf-4dir-sf08-SinRand-5s'
% '240313_SH057__no_cam_SH\run00_160812_5tf-4dir-sf08-SinRand-5s'

% EFFECT DIFFERENT SFs
% '240313_SH055__no_cam_SH\run01_140753_5sf-4dir-tf4-SinRand5s'
% '240313_SH056__no_cam_SH\run01_152943_5sf-4dir-tf4-SinRand5s'
% '240313_SH057__no_cam_SH\run01_164309_5sf-4dir-tf4-SinRand5s'
clear all
newdirs={'240311_SH057__no_cam_SH','240311_SH056__no_cam_SH'}
% sessions={'run00_170314_12dir-sf08-tf4SinSeq5s', 'run00_154140_12dir-sf08-tf4SinSeq5s'}
% sessions={'run01_172413_12dir-sf08-tf4SinSeq4s', 'run01_160616_12dir-sf08-tf4SinSeq4s'}
% sessions={'run02_174356_12dir-sf08-tf4SinSeq3s', 'run02_162601_12dir-sf08-tf4SinSeq3s'}
sessions={'run03_180209_12dir-sf08-tf4SinSeq2s', 'run03_164335_12dir-sf08-tf4SinSeq2s'}

i=2
%% read mat file with pupil trace and timestamps

path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'

for iAn=1:size(newdirs,2)
    
    newdir=newdirs{iAn};
    session=sessions{iAn};
    filepathanalysis=fullfile(path_data,  newdir, session);

    expt_name = [filepathanalysis, '\pupil_epochs.mat'];
    clear pupil_epochs;
    load(expt_name);

    pupil_trials_stim{iAn}=pupil_epochs.pupil_trials_stim;
    pupil_av_stim{iAn}=nanmean(pupil_epochs.pupil_trials_stim,2);
    pupil_trials_entire_epochs{iAn}=pupil_epochs.pupil_trials_entire_epochs;
    pupil_av_entire_epochs{iAn}=nanmean(pupil_epochs.pupil_trials_entire_epochs,2);

    pupil_01=squeeze(nanmean(pupil_av_stim{iAn}(:,:,11:22),3));
    pupil_02=squeeze(nanmean(pupil_av_stim{iAn}(:,:,end-11),3));
    delta_pupil{iAn}=(pupil_02-pupil_01)./pupil_01;

end

delta_pupil_array=cell2mat(delta_pupil);
nasal_dirs=[1,2,3,4,11,12];
temporal_dirs=[5,6,7,8,9,10];
[p,h]=ranksum(delta_pupil_array(nasal_dirs), delta_pupil_array(temporal_dirs));
% 0.6991 5s
% 0.9372 4s
%%
% tmp=cell2mat(pupil_data_av)
entire_epoch_array=cell2mat(pupil_av_entire_epochs);

stim_directions=[0:30:330];

% 0 - 5s ; 1 - 4s; 2 - 3s; 3 - 2s
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
%
clf
nAn=size(entire_epoch_array,2);


for iStim=1:12
ax_avpupil_stims(iStim) = axes('position',[.045+.08*(iStim-1),.2,.05,.1],'units','normalized');
end

for istim=1:12
axes(ax_avpupil_stims(istim))

pupil_animals_cell=squeeze(entire_epoch_array(istim,:,:))
%
% Calculate the mean and standard deviation across columns
meanData = nanmean(pupil_animals_cell);
stdData = nanstd(pupil_animals_cell);

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(nAn);
semData = stdData ;

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);

% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on

plot([0,0],[0,1],'-k','linewidth',3)

hold on;
plot([110,110],[-2,2],'--k')
hold on;
plot([110+i*22,110+i*22],[-2,2],'--k') % 5s
axis tight
ylim([-2,2])
ylabel('Pupil zscored', 'FontSize',16,'Color','k');
set(gca,'xtick',[],'ytick',[],'YColor', 'none','XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
xlim([44, length(meanData)-44])
title(sprintf('%.f\n%', stim_directions(istim)))

end


set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[session,'_averaged_epochs']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);











