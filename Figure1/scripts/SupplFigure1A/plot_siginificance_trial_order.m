% statistical test for single trial
% plot heatmaps with statistical significance

clear all
% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
session_intrest=1 % 1 - first session of the animal with stimulus

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

for itrial=1:15
    clear index_trials
    index_trials=intersect(find(new_pupil_data.raw_relative_diam_delta(:,3)==itrial), find(new_pupil_data.raw_relative_diam_delta(:,2)==session_intrest))

    selected_data{itrial}=new_pupil_data.raw_relative_diam_delta(index_trials,:);
    
    clear for_test_nasal_data
    clear for_test_temporal_data
    
    for_test_nasal_data=new_pupil_data.raw_relative_diam_delta(index_trials,3+order_nasal);
    for_test_temporal_data=new_pupil_data.raw_relative_diam_delta(index_trials,3+order_temporal);

    [pval(itrial), h0(itrial)]=signrank(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');

end

%
for iAn=1:length(pval)
    
    if pval(iAn)>0.05
        level_sign{iAn}='ns';
    elseif pval(iAn)<=0.05 &  pval(iAn)>0.01
        level_sign{iAn}='*';
    elseif pval(iAn)<=0.01 &  pval(iAn)>0.001    
        level_sign{iAn}='**';
    elseif pval(iAn)<=0.001 &  pval(iAn)>0.0001    
        level_sign{iAn}='***';
    elseif pval(iAn)<=0.0001    
        level_sign{iAn}='****';
    end

end

%%


%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS
clear for_test_nasal_data
clear for_test_temporal_data

itrial=1

mean_pupil_delta=100*selected_data{itrial}(:,4:end);
for_test_nasal_data=mean_pupil_delta(:,order_nasal);
for_test_temporal_data=mean_pupil_delta(:,order_temporal);
[pval_still, h0_still, stats_still]=signrank(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');

% Calculate the mean and standard deviation across columns
meanData = nanmean(mean_pupil_delta);
stdData = nanstd(mean_pupil_delta);

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(size(mean_pupil_delta, 1));

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);

figure(1)
% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil size change (%)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[-10, 30.1]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals_Trial01.pdf']);

%%
for itrial=1:10
mean_pupil_delta=100*selected_data{itrial}(:,4:end);
meanData_imagesc(itrial,:) = nanmean(mean_pupil_delta);
end

%%

fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);

imagesc(meanData_imagesc)
axis tight
axis square
cmap=colormap(redblue)
set(gca,'CLim',[-50 50])
ylabel('#Trials', 'FontSize',16)
xlabel('Directions (deg)','FontSize',16)
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330])
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
ylabels_animal_id={'#1','#2','#3','#4','#5','#6','#7','#8','#9','#10','#11','#12','#13'};
ylabels_animal_id=level_sign;
yticks_position=[0.5:1:13]+0.5;
set(gca,'ytick',yticks_position,'yticklabels',ylabels_animal_id);

cb=colorbar;
cb.Position = [.1,.082,.015,.05] 
caxis([-20 20]);

set(cb,'tickdir','out','fontsize',14,'ticklength',get(cb,'ticklength')*4);

box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_raw_relative_delta_pupil_1stsession_#trials.pdf']);

%%
newdir=new_pupil_data.newdir;
[animal_id_list, animal_id_index]=sort(new_pupil_data.animal_id);

newdir_first_session={'170729_BV201_2P_JC\run01_ori12_L4',
'170729_JC027_2P_JC\run02_ori12_V1'
'140623_KS092_2P_KS\run01_ori_ds_V1'
'140623_KS093_2P_KS\run01_ori_ds_V1'
'140810_KS103_2P_KS\run01_ori_ds_V1_full'
'150915_KS145_2P_KS\run01_ori_ds_V1_full'
'160630_KS164_2P_KS\run03_ori12_V1_awake'
'160916_KS166_2P_KS\run01_behav_not_rand'
'160630_KS167_2P_KS\run03_ori12_V1_awake'
'160926_KS168_2P_KS\run01_behav_not_rand'
'160926_KS169_2P_KS\run01_behav_not_rand'
'170110_KS173_2P_KS\run03_ori12_V1_awake'
'170110_KS174_2P_KS\run03_ori12_V1_awake'};

%%

for iSess=1:length(newdir_first_session);

    index=find(strcmp(newdir, newdir_first_session{iSess})==1); 
    tmp_data=new_pupil_data.trials_raw_diam_diff_relative_stimulation{index};
    ntrials=size(tmp_data,1)
    clear array_trial
    array_trial=nan(ntrials,15); 
    
    first_session_index(iSess)=tmp_data

    array_trial(1:ntrials,1)=new_pupil_data.animal_id(index);
    array_trial(1:ntrials,2)=1;
    array_trial(1:ntrials,3)=1:1:ntrials;
    array_trial(1:ntrials,4:end)=tmp_data;
    
    appendedArray = [array1; array2];
end


%%


% Calculate the average pupil size for each trial and animal
averagePupilSize = grpstats(data, {'animalID', 'trialNumber'}, 'mean', 'DataVars', 'pupilSize');

% Group the data by trial and calculate the mean pupil size
meanPupilSizeByTrial = grpstats(averagePupilSize, 'trialNumber', 'mean', 'DataVars', 'mean_pupilSize');

% Create a matrix to store the average pupil size for each trial
averagePupilMatrix = reshape(meanPupilSizeByTrial.mean_pupilSize, [], 12);

% Perform linear regression for each direction separately
slopeCoefficients = zeros(1, 12);
pValues = zeros(1, 12);

for i = 1:12
    model = fitlm(meanPupilSizeByTrial, ['mean_pupilSize_' num2str(i) ' ~ trialNumber']);
    slopeCoefficients(i) = model.Coefficients.Estimate(2);
    pValues(i) = model.Coefficients.pValue(2);
end

% Display the results
disp('Slope coefficients:');
disp(slopeCoefficients);

disp('p-values:');
disp(pValues);

%%

% Assuming your data is in a table named 'data' with columns 'pupilSize', 'animalID', and 'trialNumber'

% Calculate the average pupil size for each trial and animal
averagePupilSize = grpstats(data, {'animalID', 'trialNumber'}, 'mean', 'DataVars', 'pupilSize');

% Group the data by trial and calculate the mean pupil size
meanPupilSizeByTrial = grpstats(averagePupilSize, 'trialNumber', 'mean', 'DataVars', 'mean_pupilSize');

% Perform linear regression
model = fitlm(meanPupilSizeByTrial, 'mean_pupilSize ~ trialNumber');

% Display the regression results
disp(model);

% Assess the significance of the slope coefficient
disp('Slope coefficient p-value:');
disp(model.Coefficients.pValue(2));

%%
size([new_pupil_data.trials_raw_diam_diff_relative_stimulation{first_session_index}])

