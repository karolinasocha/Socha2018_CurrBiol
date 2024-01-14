%%
% statistical test for single trial
% plot heatmaps with statistical significance

clear all
% this load processed data for checking statistical tests
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_0_05.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_05_1.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_1_15.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_15_2.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_2_25.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_25_3.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_3_35.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data_time_35_4.mat')
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')

% pvals for individual animals
clear h0
clear pval
clear ci
clear stats

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

animal_id_list=unique(new_pupil_data.animal_id)
for iAn=1:length(animal_id_list)

    indecies=find(new_pupil_data.animal_id==animal_id_list(iAn))
    appendedArray=[];
    
    for ii=1:length(indecies)
        clear tmp
        %tmp=new_pupil_data.trials_median_diam_stims_reordered_offset{indecies(ii)};
        tmp=new_pupil_data.trials_raw_diam_diff_relative_stimulation{indecies(ii)};
        
        appendedArray = vertcat(appendedArray, tmp);
    end
    
    data_nasal=appendedArray(:,order_nasal);
    data_temporal=appendedArray(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    %[pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
     %[pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    [pval(iAn), h0(iAn), stats{iAn}]=ranksum(data_nasal(:), data_temporal(:),'Tail','right');
    pupil_animals_cell{iAn}=squeeze(nanmean(appendedArray,1))
    
end

number_increased_diam_animals=length(find(pval<0.05));

% scatter plot COLOR CODED ANIMAL

colors=jet(13);

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

av_delta_animals=cell2mat(pupil_animals_cell')
figure
for iAn=1:length(pupil_animals_cell)
    clear data_temporal
    clear data_nasal
    
    data_nasal=pupil_animals_cell{iAn}(order_nasal)
    data_temporal=pupil_animals_cell{iAn}(order_temporal)
    hold on
%   plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', colors(iAn,:),'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)

    data_nasal_all{iAn}=nanmean(data_nasal);
    data_temporal_all{iAn}=nanmean(data_temporal);
end

hold on
plot([-100 100], [-100 100],'--k');

ylim([-10, 70])
xlim([-10, 40])

for_test_nasal_data=cell2mat(data_nasal_all');
for_test_temporal_data=cell2mat(data_temporal_all');

% for plots calculate average point 
mean_all_nasal=100*nanmean(for_test_nasal_data(:),1);
sem_all_nasal=100*nanstd(for_test_nasal_data(:),1)./sqrt(size(for_test_nasal_data, 1));

mean_all_temporal=100*nanmean(for_test_temporal_data(:),1);
sem_all_temporal=100*nanstd(for_test_temporal_data(:),1)./sqrt(size(for_test_temporal_data, 1));

%[pval, h0, stats]=signrank(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');
[pval, h0, stats]=ranksum(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');


% pval= 0.0324 1-1.5
% pval= 0.0045 1.5-2
% pval =0.0017 2-2.5
% pval= 8.7939e-04 2.5-3
% pval= 5.1527e-04 3-3.5
% pval=2.9530e-04 3.5-4
% pval= 7.3860e-05 4.5-5
%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS
% meanData_25_3;
% meanData_3_35;

mean_pupil_delta=100*cell2mat(pupil_animals_cell');

% Calculate the mean and standard deviation across columns
meanData = nanmean(mean_pupil_delta);
stdData = nanstd(mean_pupil_delta);

save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_05_1.mat','meanData')

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(size(mean_pupil_delta, 1));
save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_05_1.mat','semData')

% Define x-axis values (e.g., assuming 5 data points)
x = 1:numel(meanData);
figure(1)
clf
% Plot the mean data
plot(x, meanData, 'k', 'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil size change (%)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[-5, 30.1]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
% filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\'

print(gcf,'-dpdf',[filepathanalysis, 'Figure1E_average_raw_diam_relative_delta_pupil_SEM_animals_05_1.pdf']);

% calculate mean + SEM
% meanData =[16.1359   12.2291    8.1056    4.2092   -0.3538    1.4475    0.2671   -1.1895   -0.1307    3.9853    4.9183   23.9914]
%semData=[4.1998    2.1429    2.2035    1.7687    1.6019    3.4430    1.3093    1.1727    1.0122    1.8115    1.7518    5.5431 ]

% 0deg: 16.1359 +- 4.1998 (across 13 animals)
% 90deg: 4.2092 +- 1.7687 (across 13 animals)
% 180deg: 0.2671 +- 1.3093 (across 13 animals)
% 270deg: 3.9853 +- 1.8115 (across 13 animals)

%%
