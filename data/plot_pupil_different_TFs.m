% '240313_SH055__no_cam_SH\run00_133318_5tf-4dir-sf08-SinRand-5s' 
% '240313_SH056__no_cam_SH\run00_145132_5tf-4dir-sf08-SinRand-5s'
% '240313_SH057__no_cam_SH\run00_160812_5tf-4dir-sf08-SinRand-5s'
%% read mat file with pupil trace and timestamps

clear all

newdirs={'240313_SH055__no_cam_SH','240313_SH056__no_cam_SH','240313_SH057__no_cam_SH'}
sessions={'run00_133318_5tf-4dir-sf08-SinRand-5s', 'run00_145132_5tf-4dir-sf08-SinRand-5s','run00_160812_5tf-4dir-sf08-SinRand-5s'}

path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'


for iAn=1:size(newdirs,2)
    newdir=newdirs{iAn};
    session=sessions{iAn};
    filepathanalysis=fullfile(path_data,  newdir, session);

    expt_name = [filepathanalysis, '\pupil_epochs.mat'];
    clear pupil_epochs;
    load(expt_name);

    pupil_trials_stim{iAn}=pupil_epochs.pupil_trials_stim;
    pupil_av_stim{iAn}=nanmean(pupil_epochs.pupil_trials_stim,3);
    pupil_trials_entire_epochs{iAn}=pupil_epochs.pupil_trials_entire_epochs;
    pupil_av_entire_epochs{iAn}=nanmean(pupil_epochs.pupil_trials_entire_epochs,3);

%     pupil_01=squeeze(nanmean(pupil_av_stim{iAn}(:,:,:,11:22),4));
%     pupil_02=squeeze(nanmean(pupil_av_stim{iAn}(:,:,:,end-11),4));
%     delta_pupil{iAn}=(pupil_02-pupil_01)./pupil_01;

end
%%

pupil_av_entire_epochs_array = [];

for i = 1:length(pupil_av_entire_epochs)
    pupil_av_entire_epochs_array = cat(5, pupil_av_entire_epochs_array, pupil_av_entire_epochs{i});
end 
 
 %%
stim_directions=pupil_epochs.dirs_stimulus;
stim_tfs=pupil_epochs.tfs_stimulus;
 

 %%
% Assuming stim_directions and stim_sfs define the grid size (4x5 in this case)
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
for istim = 1:length(stim_directions) % Assuming this iterates over rows
    for itf = 1:length(stim_tfs) % Assuming this iterates over columns
        % Calculate position
        % Adjust the spacing values accordingly to fit the plots within the figure window
        left = .07 + .1 * (itf - 1); % Adjust horizontal spacing based on column
        bottom = .8 - .15 * (istim - 1); % Adjust vertical spacing based on row
        width = .05; % Fixed width
        height = .1; % Fixed height
        
        % Create axes for the plot
        ax_avpupil_nasal(istim,itf) = axes('position',[left, bottom, width, height],'units','normalized');
    end
end

 %%
i=5
nAn=size(pupil_av_entire_epochs_array,5);

 for istim=1:length(stim_directions)
    for itf = 1:length(stim_tfs)
    axes(ax_avpupil_nasal(istim,itf))

    clear pupil_animals_cell
    pupil_animals_cell=squeeze(pupil_av_entire_epochs_array(istim,itf,:,:,:));
%
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

title(sprintf('%.fdeg\n%.fHz', stim_directions(istim), stim_tfs(itf)));


end

 end
 
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[session,'_averaged_epochs_differentTFs']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);
 
 