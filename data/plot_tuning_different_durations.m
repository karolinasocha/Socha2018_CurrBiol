clear all

% EFFECT DIFFERENT DURATION
% '240311_SH056__no_cam_SH\run00_154140_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH056__no_cam_SH\run01_160616_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH056__no_cam_SH\run02_162601_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH056__no_cam_SH\run03_164335_12dir-sf08-tf4SinSeq2s' - 2s duration

% '240311_SH057__no_cam_SH\run00_170314_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH057__no_cam_SH\run01_172413_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH057__no_cam_SH\run02_174356_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH057__no_cam_SH\run03_180209_12dir-sf08-tf4SinSeq2s' - 2s duration


% 'SH0550304\run00_112038_12dir-sf08-tf-4SinRan5s'
% 'SH0550304\run01_120253_12dir-sf08-tf-4SinSeq2s'
% 'SH0550304\run02_121818_12dir-sf08-tf-4SinSeq3s'
% 'SH0550304\run03_123655_12dir-sf08-tf-4SinSeq4s'
% 'SH0550304\run04_130132_12dir-sf08-tf-4SinSeq5s'

clear all
path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'
newdirs={'240311_SH056__no_cam_SH','240311_SH057__no_cam_SH','SH0550304'}
% sessions={'run03_164335_12dir-sf08-tf4SinSeq2s','run03_180209_12dir-sf08-tf4SinSeq2s','run01_120253_12dir-sf08-tf-4SinSeq2s'}
% sessions={'run02_162601_12dir-sf08-tf4SinSeq3s', 'run02_174356_12dir-sf08-tf4SinSeq3s','run02_121818_12dir-sf08-tf-4SinSeq3s'}
% sessions={'run01_160616_12dir-sf08-tf4SinSeq4s', 'run01_172413_12dir-sf08-tf4SinSeq4s','run03_123655_12dir-sf08-tf-4SinSeq4s'}
sessions={'run00_154140_12dir-sf08-tf4SinSeq5s', 'run00_170314_12dir-sf08-tf4SinSeq5s','run00_112038_12dir-sf08-tf-4SinRan5s'}

path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'

i=5
for iAn=1:size(newdirs,2)
    newdir=newdirs{iAn};
    session=sessions{iAn};
    filepathanalysis=fullfile(path_data,  newdir, session);

    expt_name = [filepathanalysis, '\pupil_epochs.mat'];
    clear pupil_epochs;
    load(expt_name);

    pupil_trials_stim{iAn}=pupil_epochs.pupil_trials_stim;
    pupil_av_stim_end{iAn}=squeeze(nanmean(nanmean(pupil_epochs.pupil_trials_stim(:,:,end-11:end),2),3));
    pupil_av_stim_start{iAn}=squeeze(nanmean(nanmean(pupil_epochs.pupil_trials_stim(:,:,11:22),2),3));
    pupil_av_stim{iAn}=nanmean(pupil_epochs.pupil_trials_stim,3);
    pupil_trials_entire_epochs{iAn}=pupil_epochs.pupil_trials_entire_epochs;
    pupil_av_entire_epochs{iAn}=nanmean(pupil_epochs.pupil_trials_entire_epochs,2);
    delta_pupil{iAn}=(pupil_av_stim_end{iAn}-pupil_av_stim_start{iAn})
%     pupil_01=squeeze(nanmean(pupil_av_stim{iAn}(:,:,:,11:22),4));
%     pupil_02=squeeze(nanmean(pupil_av_stim{iAn}(:,:,:,end-11),4));
%     delta_pupil{iAn}=(pupil_02-pupil_01)./pupil_01;

end
%%
pupil_animals_cell = cell2mat(pupil_av_stim_end);
pupil_delta_cell = cell2mat(delta_pupil);

pupil_av_entire_epochs_array=cell2mat(pupil_av_entire_epochs);
nAn=size(pupil_av_stim_end,2);

% pupil_av_entire_epochs_array = [];
% 
% for i = 1:length(pupil_av_entire_epochs)
%     pupil_av_entire_epochs_array = cat(4, pupil_av_entire_epochs_array, pupil_av_entire_epochs{i});
% end 
%  

%
stim_directions=pupil_epochs.dirs_stimulus;
stim_tfs=pupil_epochs.tfs_stimulus;
stim_sfs=pupil_epochs.sfs_stimulus;
%
% Assuming stim_directions and stim_sfs define the grid size (4x5 in this case)
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);

for iStim=1:12
ax_epochpupil(iStim) = axes('position',[.045+.08*(iStim-1),.2,.05,.1],'units','normalized');
end

ax_avpupil = axes('position',[0.2, 0.6, 0.1, 0.15],'units','normalized');
ax_deltapupil = axes('position',[0.4, 0.6, 0.1, 0.15],'units','normalized');

axes(ax_avpupil)
% Calculate the mean and standard deviation across columns
meanData = squeeze(nanmean(pupil_animals_cell,2))';
stdData = nanstd(pupil_animals_cell');

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(nAn);
% semData = stdData ;

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);

% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
hold on
% plot([0,0],[-.5,0],'-k','linewidth',3)
ylim([-1,1])
ylabel('Pupil zscored', 'FontSize',16,'Color','k');
% set(gca,'xtick',[],'ytick',[],'YColor', 'none','XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);

box off

axes(ax_deltapupil)
% Calculate the mean and standard deviation across columns
meanData = squeeze(nanmean(pupil_delta_cell,2))';
stdData = nanstd(pupil_delta_cell');

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(nAn);
% semData = stdData ;

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);

% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
hold on
% plot([0,0],[-.5,0],'-k','linewidth',3)
ylim([-1,1])
ylabel('delta pupil', 'FontSize',16,'Color','k');
% set(gca,'xtick',[],'ytick',[],'YColor', 'none','XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
set(gca,'xtick',[1:3:12],'xticklabel',[0:90:330],'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);

box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[session,'_averaged_epochs_duration_3mice_3s']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);


%% %
for istim=1:12
axes(ax_epochpupil(istim))


% Calculate the mean and standard deviation across columns
meanData = squeeze(nanmean(pupil_av_entire_epochs_array(istim,:,:),2))';
stdData = squeeze(nanstd(pupil_av_entire_epochs_array(istim,:,:)));
% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(nAn);
% semData = stdData ;

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);

% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;
% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
hold on
plot([44,44],[0,1],'-k','linewidth',3)

hold on;
plot([110,110],[-2,2],'--k')
hold on;
plot([110+i*22,110+i*22],[-2,2],'--k') % 5s
axis tight
ylim([-1,1])
ylabel('Pupil zscored', 'FontSize',16,'Color','k');
set(gca,'xtick',[],'ytick',[],'YColor', 'none','XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
xlim([44, length(meanData)-44])
title(sprintf('%.f\n%', stim_directions(istim)))
end



set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[session,'_averaged_epochs_duration_3mice_2']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);







