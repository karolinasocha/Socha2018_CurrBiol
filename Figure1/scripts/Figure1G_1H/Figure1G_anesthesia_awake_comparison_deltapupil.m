%% code written 21-07-2023 by KSocha
% correct version to calculate mean +- SEM for anesthesia vs awake
% Figure 1G
clear all

newdir_diameter_anesthesia={'181108_KS600_ephys_behav_JC\run01_direction_12dir_4Hz_ansth',...
    '181108_KS602_ephys_behav_JC\run01_direction_12dir_4Hz_anesth',...
    '170110_KS173_2P_KS\run03_ori12_V1_anesthesia'};

newdir_diameter_awake={'181108_KS600_ephys_behav_JC\run02_direction_12dir_4Hz',...
     '181108_KS602_ephys_behav_JC\run02_direction_12dir_4Hz',...
     '170110_KS173_2P_KS\run03_ori12_V1_awake_after'};
 
 %%
 
 
 %% load ANESTHESIA data
 
for iAn=1:length(newdir_diameter_anesthesia)
    clear diameter_data_recalculated_revision
    filepathanalysis=['J:\data\analysis\' newdir_diameter_anesthesia{iAn}]; 
    load([filepathanalysis, '\diameter_data_recalculated_revision.mat']);

    diameter_data_recalculated_anesthesia{iAn}=diameter_data_recalculated_revision;
    try
    ddat_interp_area_anesthesia{iAn}=diameter_data_recalculated_anesthesia{iAn}.diameter_area.ddat_interp_area;
    end
    try
        ddat_interp_area_anesthesia{iAn}=diameter_data_recalculated_anesthesia{iAn}.diameter_area.ddat_interp_area_epochs
    end
    pre_frames_interpolated=diameter_data_recalculated_anesthesia{1}.pre_frames_interpolated;
    stim_frames_inerpolated=diameter_data_recalculated_anesthesia{1}.stim_frames_inerpolated;
    
    frameRate_interpolated=diameter_data_recalculated_anesthesia{1}.frameRate_interpolated;
    onset=round(frameRate_interpolated*0.5);
    offset=round(frameRate_interpolated*1);
    
    ddat_interp_area_stim_anesthesia{iAn}=ddat_interp_area_anesthesia{iAn}(pre_frames_interpolated:stim_frames_inerpolated,:,:);

    ddat_interp_area_stim_onset_anesthesia{iAn}=squeeze(nanmean(ddat_interp_area_stim_anesthesia{iAn}(onset:offset,:,:),1));
    ddat_interp_area_stim_offset_anesthesia{iAn}=squeeze(nanmean(ddat_interp_area_stim_anesthesia{iAn}(end-onset:end,:,:),1));
    
    delta_pupil_stim_anesthesia{iAn}=(ddat_interp_area_stim_offset_anesthesia{iAn}-ddat_interp_area_stim_onset_anesthesia{iAn})./ddat_interp_area_stim_onset_anesthesia{iAn}
   
    av_anasthesia_raw_diam_diff_relative_stimulation{iAn}=nanmean(delta_pupil_stim_anesthesia{iAn},1) ;


end

%% load AWAKE
 
for iAn=1:length(newdir_diameter_awake)
    clear diameter_data_recalculated_revision
    filepathanalysis=['J:\data\analysis\' newdir_diameter_awake{iAn}]; 
    load([filepathanalysis, '\diameter_data_recalculated_revision.mat']);

    diameter_data_recalculated_awake{iAn}=diameter_data_recalculated_revision;
    try
    ddat_interp_area_awake{iAn}=diameter_data_recalculated_awake{iAn}.diameter_area.ddat_interp_area;
    end
    try
        ddat_interp_area_awake{iAn}=diameter_data_recalculated_awake{iAn}.diameter_area.ddat_interp_area_epochs
    end
    pre_frames_interpolated=diameter_data_recalculated_awake{1}.pre_frames_interpolated;
    stim_frames_inerpolated=diameter_data_recalculated_awake{1}.stim_frames_inerpolated;
    
    frameRate_interpolated=diameter_data_recalculated_awake{1}.frameRate_interpolated;
    onset=round(frameRate_interpolated*0.5);
    offset=round(frameRate_interpolated*1);
    
    ddat_interp_area_stim_awake{iAn}=ddat_interp_area_awake{iAn}(pre_frames_interpolated:stim_frames_inerpolated,:,:);

    ddat_interp_area_stim_onset_awake{iAn}=squeeze(nanmean(ddat_interp_area_stim_awake{iAn}(onset:offset,:,:),1));
    ddat_interp_area_stim_offset_awake{iAn}=squeeze(nanmean(ddat_interp_area_stim_awake{iAn}(end-onset:end,:,:),1));
    
    delta_pupil_stim_awake{iAn}=(ddat_interp_area_stim_offset_awake{iAn}-ddat_interp_area_stim_onset_awake{iAn})./ddat_interp_area_stim_onset_awake{iAn}
   
    av_awake_raw_diam_diff_relative_stimulation{iAn}=nanmean(delta_pupil_stim_awake{iAn},1) ;

end

%%

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

mean_pupil_delta_anesthesia=100*cell2mat(av_anasthesia_raw_diam_diff_relative_stimulation');
mean_pupil_delta_anesthesia=mean_pupil_delta_anesthesia(:,1:12);
% Calculate the mean and standard deviation across columns
meanData_ansth = nanmean(mean_pupil_delta_anesthesia);
stdData_ansth = nanstd(mean_pupil_delta_anesthesia);

% Calculate the Standard Error of the Mean (SEM)
semData_ansth = stdData_ansth ./ sqrt(size(mean_pupil_delta_anesthesia, 1));

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData_ansth);
figure(1)
% Plot the mean data
plot(x, meanData_ansth, 'b', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_ansth - semData_ansth, fliplr(meanData_ansth + semData_ansth)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% AWAKE

mean_pupil_delta_awake=100*cell2mat(av_awake_raw_diam_diff_relative_stimulation');
mean_pupil_delta_awake=mean_pupil_delta_awake(:,1:12);

% Calculate the mean and standard deviation across columns
meanData_awake = nanmean(mean_pupil_delta_awake);
stdData_awake = nanstd(mean_pupil_delta_awake);

% Calculate the Standard Error of the Mean (SEM)
semData_awake = stdData_awake ./ sqrt(size(mean_pupil_delta_awake, 1));

plot(x, meanData_awake, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_awake - semData_awake, fliplr(meanData_awake + semData_awake)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil size change (%)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[-5, 40]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

%print(gcf,'-dpdf',[filepathanalysis, 'Figure1G_3mice_average_anesthesia_awake_relative_delta_pupil_SEM_animals.pdf']);
%
%% PVAL - TEST

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

data_nasal_awake=mean_pupil_delta_awake(:,order_nasal);
data_temporal_awake=mean_pupil_delta_awake(:,order_temporal);
%[pval_awake, h0_awake, stats_awake]=signrank(data_nasal_awake(:), data_temporal_awake(:),'Tail','right');
[pval_awake, h0_awake, stats_awake]=ranksum(data_nasal_awake(:), data_temporal_awake(:),'Tail','right');

pval_awake % pval=0.0118

data_nasal_anesthesia=mean_pupil_delta_anesthesia(:,order_nasal);
data_temporal_anesthesia=mean_pupil_delta_anesthesia(:,order_temporal);
[pval_anesthesia, h0_anesthesia, stats_anesthesia]=signrank(data_nasal_anesthesia(:), data_temporal_anesthesia(:),'Tail','right');
[pval_anesthesia, h0_anesthesia, stats_anesthesia]=ranksum(data_nasal_anesthesia(:), data_temporal_anesthesia(:),'Tail','right');

pval_anesthesia % pval=0.2040