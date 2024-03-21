

clear all

path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'
newdir='240314_VO050'
session='run00_female'
filepathanalysis=fullfile(path_data,  newdir, session); 
expt_name = [filepathanalysis, '\KSdata.mat'];

experiment_data=load(expt_name);
fps=experiment_data.fsl
fps=1
stimuli_params=experiment_data.stimuli
pupil_data=experiment_data.timepar
pupil_timestamps=pupil_data(:,1);
pupil_area=pupil_data(:,2);
pupil_area_zscore = bsxfun(@rdivide, bsxfun(@minus, pupil_area, nanmean(pupil_area)) , nanstd(pupil_area));
pupil_area_zscore2=zscore(pupil_area);
stimulusOnset=experiment_data.vstim_onoff(:,1);
stimulusOffset=experiment_data.vstim_onoff(:,2);
%%

figure, 
subplot(3,1,1)
plot(pupil_area,'k');
hold on
set(gca,'xtick',[],'XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
axis tight
ylim([0,4000])
title(sprintf('Raw\n%s', newdir))

subplot(3,1,2)
smoothedPupil = smoothdata(pupil_area,'sgolay', 10);
plot(smoothedPupil,'k');
hold on
set(gca,'xtick',[],'XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
axis tight
ylim([0,4000])
title(sprintf('SmoothRaw\n%s', newdir))

subplot(3,1,3)
plot(pupil_area_zscore,'k');
hold on
set(gca,'xtick',[],'XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
axis tight
ylim([-4,4])
title(sprintf('zscored\n%s', session))

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[newdir,'_',session,'_rawtrace']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);
%%
stim_directions=unique(stimuli_params(:,1));
stim_tfs=unique(stimuli_params(:,3));
stim_sfs=unique(stimuli_params(:,2));
ntrials= length(find(stimuli_params(:,1)==0));

stim_duration=5;
prestim_duration=5;
poststim_duration=5;

pupil_trials_stim=nan(length(stim_directions), ntrials,22*stim_duration); % 110 for 5s
pupil_trials_poststim=nan(length(stim_directions), ntrials,22*stim_duration);
pupil_trials_epochs=nan(length(stim_directions), ntrials,110+22*stim_duration); % 220 for 5s
pupil_trials_entire_epochs=nan(length(stim_directions), ntrials,220+22*stim_duration); % 330 for 5s

for istim=1:length(stim_directions)
        
    indx=find(stimuli_params(:,1)==stim_directions(istim));
    stimOnset_indx=stimulusOnset(indx);
    stimOffset_indx=stimulusOffset(indx) 
    stimOffset_indx_fixed=5+stimOnset_indx; 
    
    for itrial=1:length(stimOnset_indx)
        
        target_timestamp = stimOnset_indx(itrial);
        differences = abs(pupil_timestamps - target_timestamp);
        [~, nearestIndex_onset] = min(differences);
        
        target_timestampOffset = stimOffset_indx(itrial);
        differences2 = abs(pupil_timestamps - target_timestampOffset);
        [~, nearestIndex_offset] = min(differences2); 
        
        target_timestamppre = stimOnset_indx(itrial)-prestim_duration*fps;
        differences3 = abs(pupil_timestamps - target_timestamppre);
        [~, nearestIndex_onsetPre] = min(differences3);
        
        target_timestamppost = stimOffset_indx(itrial)+poststim_duration*fps;
        differences4 = abs(pupil_timestamps - target_timestamppost);
        [~, nearestIndex_offsetPost] = min(differences4);
        
        clear tmp
        tmp=pupil_area_zscore(nearestIndex_onset:nearestIndex_offset);
        pupil_trials_stim(istim, itrial,1:length(tmp))=tmp;
        
        clear tmp2
        tmp2=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_onset);
        pupil_trials_blanks(istim, itrial,1:length(tmp2))=tmp2;
        
        clear tmp3
        tmp3=pupil_area_zscore(nearestIndex_offset:nearestIndex_offsetPost);
        pupil_trials_poststim(istim, itrial,1:length(tmp3))=tmp3;
        
        clear tmp4
        tmp4=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_offset);
        pupil_trials_epochs(istim, itrial,1:length(tmp4))=tmp4;
        
        clear tmp5
        tmp5=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_offsetPost);
        pupil_trials_entire_epochs(istim, itrial,1:length(tmp5))=tmp5;

    end  
    
end

pupil_epochs.pupil_trials_entire_epochs=pupil_trials_entire_epochs
pupil_epochs.pupil_trials_blanks=pupil_trials_blanks
pupil_epochs.pupil_trials_stim=pupil_trials_stim

pupil_epochs.pupil_trials_poststim=pupil_trials_poststim
pupil_epochs.pupil_trials_epochs=pupil_trials_epochs

pupil_epochs.stimulus=stim_directions
pupil_epochs.dirs_stimulus=stim_directions
pupil_epochs.tfs_stimulus=stim_tfs
pupil_epochs.sfs_stimulus=stim_sfs
pupil_epochs.ntrials=ntrials
pupil_epochs.prestim_duration=prestim_duration
pupil_epochs.poststim_duration=poststim_duration
pupil_epochs.stim_duration=stim_duration
pupil_epochs.type='fixed sequence'
pupil_epochs.stimulus_duration_in_sec='5s'

pupil_mat_name = [filepathanalysis, '\pupil_epochs.mat'];
save(pupil_mat_name, 'pupil_epochs');

%%
% tmp=cell2mat(pupil_data_av)
entire_epoch_array=pupil_trials_entire_epochs;

stim_directions=[0:30:330];

% 0 - 5s ; 1 - 4s; 2 - 3s; 3 - 2s
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[20,20],'paperposition',[0,0,20,20]);
%
clf
ntrials=size(entire_epoch_array,2);
i=5
n_frames=30.4000 % 22

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
semData = stdData ./ sqrt(ntrials);
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
plot([n_frames*5,n_frames*5],[-2,2],'--k')
hold on;
plot([n_frames*5+i*n_frames,n_frames*5+i*n_frames],[-2,2],'--k') % 5s
hold on;
plot([n_frames*2,n_frames*2],[0,1],'k','linewidth',3)
% plot([110+i*n_frames,110+i*n_frames],[-2,2],'--k') % 5s
axis tight
ylim([-2,2])
ylabel('Pupil zscored', 'FontSize',16,'Color','k');
set(gca,'xtick',[],'ytick',[],'YColor', 'none','XColor', 'none','tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off
xlim([n_frames*2, length(meanData)-n_frames*2])
title(sprintf('%.f\n%', stim_directions(istim)))

end


set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
savefigurepath='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\figures\'
% filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\'
filename=[session,'_averaged_epochs']
print(gcf,'-dpdf',[savefigurepath, filename,'.pdf']);
print(gcf,'-dpdf',[savefigurepath, filename,'.png']);















