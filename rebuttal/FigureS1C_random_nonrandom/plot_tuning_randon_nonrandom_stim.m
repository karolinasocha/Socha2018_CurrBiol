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
diam_interpolated_av{iAn}=nanmean(diam_interpolated{iAn},2);

end

diam_epochs=cell2mat(diam_interpolated);
diam_epochs_random=cell2mat(diam_interpolated);
diam_epochs_random=cell2mat(diam_interpolated);
velo_epochs=cell2mat(velo_animals);

diam_stim_random_array=cell2mat(diam_interpolated_av);
animal01_R=nanmean(diam_stim_random_array(:,1:2,:),2);
diam_animals_R=[animal01_R, diam_stim_random_array(:,3:4,:)];

diam_stim_random_av=squeeze(nanmean(diam_animals_R(195:212,:,:),1));
diam_stim_random_av=diam_stim_random_av(:,1:12)

% diam_stim_random_av=squeeze(nanmean(diam_stim_random_array(195:212,:,:),1));
% diam_stim_random_av=diam_stim_random_av(:,1:12)

%%
newdir_nonrandom={'180823_KS501_2P_JC\run00_direction_12dir_4Hz',...
 '180824_KS501_2P_JC\run00_direction_12dir_4Hz',...
 '180823_KS503_2P_JC\run00_direction_12dir_4Hz',...
 '180824_KS503_2P_JC\run00_direction_12dir_4Hz'};

analysis_path='J:\data\analysis\';

for iAn=1:size(newdir_random,2)
   clear tmp_behav 
tmp_behav=load([analysis_path newdir_nonrandom{iAn} '\behav_data.mat']);
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
clear diam_interpolated_av

for iAn=1:size(newdir_random,2)
    clear dtime
    clear p2time
    dtime = linspace(0,1,size(diam_animals{iAn},1));
    p2time =  linspace(0,1,size(diam_animals{3},1));
for iStim = 1:13
    for iTrial=1:10
        clear tempo
        tempo=diam_animals{iAn}(:,iTrial,iStim);
        diam_interpolated{iAn}(:,iTrial,iStim)=interp1(dtime,tempo,p2time,'linear','extrap');
    end
end
diam_interpolated_av{iAn}=nanmean(diam_interpolated{iAn},2);

end


diam_epochs_non_random=cell2mat(diam_interpolated);
velo_epochs=cell2mat(velo_animals);

diam_stim_non_random_array=cell2mat(diam_interpolated_av);
animal01=nanmean(diam_stim_non_random_array(:,1:2,:),2);
diam_animals=[animal01, diam_stim_non_random_array(:,3:4,:)];

diam_stim_non_random_av=squeeze(nanmean(diam_animals(195:212,:,:),1));
diam_stim_non_random_av=diam_stim_non_random_av(:,1:12);

% diam_stim_non_random_av=squeeze(nanmean(diam_stim_non_random_array(30:60,:,:),1));
% diam_stim_non_random_av=diam_stim_non_random_av(:,1:12);


%%
%%
meanData = nanmean(diam_stim_random_av);
stdData = nanstd(diam_stim_random_av);
semData = stdData ./ sqrt(size(diam_stim_random_av, 1));

meanData_nR = nanmean(diam_stim_non_random_av);
stdData_nR = nanstd(diam_stim_non_random_av);

% Calculate the Standard Error of the Mean (SEM)
semData_nR = stdData_nR ./ sqrt(size(diam_stim_non_random_av, 1));

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);
figure(1)
clf
% Plot the mean data
plot(x, meanData, 'r', 'LineWidth', 2);
hold on;
plot(x, meanData_nR, 'k', 'LineWidth', 2);
hold on;
% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
fill([x, fliplr(x)], [meanData_nR - semData_nR, fliplr(meanData_nR + semData_nR)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil area (mm2)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[0.6, 1.2]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);

box off
title('stimulation')
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
% filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

%%
set(gcf,'paperunits','centimeters','papersize' ,[22,22],'color','w','paperposition',[0,0,21,21],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\FigureS1C_random_nonrandom\']; 
print(gcf,'-dpdf',[filepathanalysis, '\pupil_stimulus_results_3mice_nonrand_random_3mice_STIM.pdf']);
print(gcf,'-dpng',[filepathanalysis, '\pupil_stimulus_results_3mice_nonrand_random_3mice_STIM.png']);
