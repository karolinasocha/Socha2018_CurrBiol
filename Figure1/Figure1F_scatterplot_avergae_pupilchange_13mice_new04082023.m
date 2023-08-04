%%
% statistical test for single trial
% plot heatmaps with statistical significance

clear all
% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')

%% pvals for individual animals
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
    [pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    pupil_animals_cell{iAn}=squeeze(nanmean(appendedArray,1));
    
    trials_nasal{iAn}=appendedArray(:,order_nasal);
    trials_temporal{iAn}=appendedArray(:,order_temporal);
    
end

number_increased_diam_animals=length(find(pval<0.05));
number_increased_diam_animals
%%

for iAn=1:length(trials_nasal)
    
    trials_nasal_mean{iAn}=nanmean(trials_nasal{iAn});
    trials_temporal_mean{iAn}=nanmean(trials_temporal{iAn});
    num_trials = size(trials_nasal_mean{iAn}, 2);
    alpha = 0.05; % 95% confidence level

    % Calculate mean and standard error of the mean (SEM) for each condition
    mean_values_nasal{iAn} = nanmean(trials_nasal_mean{iAn});
    sem_values_nasal{iAn} = std(trials_nasal_mean{iAn}) / sqrt(num_trials);

    mean_values_temporal{iAn} = nanmean(trials_temporal_mean{iAn});
    sem_values_temporal{iAn} = std(trials_temporal_mean{iAn}) / sqrt(num_trials);
    
    t_value = tinv(1 - alpha/2, num_trials - 1);

    % Calculate confidence intervals
    lower_ci_nasal{iAn} = mean_values_nasal{iAn} - t_value * sem_values_nasal{iAn};
    upper_ci_nasal{iAn} = mean_values_nasal{iAn} + t_value * sem_values_nasal{iAn};

    % Calculate confidence intervals
    lower_ci_temporal{iAn} = mean_values_temporal{iAn} - t_value * sem_values_temporal{iAn};
    upper_ci_temporal{iAn} = mean_values_temporal{iAn} + t_value * sem_values_temporal{iAn};
      
end
%% WITH SEM
colors=jet(13);
colors_transparency=colors;
colors_transparency=colors_transparency;
figure(2)
clf
for iAn=1:length(trials_temporal_mean)
    data_temporal=trials_temporal_mean{iAn};
    data_nasal=trials_nasal_mean{iAn};
    hold on
    %plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    
    %plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', colors_transparency(iAn,:),'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    s1=scatter(100*data_temporal(:),100*data_nasal(:),'o','MarkerFaceColor',colors(iAn,:),'MarkerEdgeColor','k');
    s1.MarkerFaceAlpha = .5;
    
    % Create a plot with vertical error bars and transparency (alpha = 0.5)
    %x = 1:size(data_example, 2);
    errorbar(100*mean_values_temporal{iAn}, 100*(mean_values_nasal{iAn} - sem_values_nasal{iAn}),100*(mean_values_nasal{iAn} + sem_values_nasal{iAn}),100*mean_values_nasal{iAn}, 100*(mean_values_temporal{iAn} - sem_values_temporal{iAn}),100*( mean_values_temporal{iAn} + sem_values_temporal{iAn}), 'o', 'MarkerFaceColor', colors(iAn,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5,'Color','k', 'MarkerSize', 12);
    hold on;

    % Create a plot with horizontal  and vertical CI
    %errorbar(100*mean_values_temporal{iAn}, 100*mean_values_nasal{iAn},100*(mean_values_nasal{iAn} - lower_ci_nasal{iAn}),100*( upper_ci_nasal{iAn} - mean_values_nasal{iAn}), 100*(mean_values_temporal{iAn} - lower_ci_temporal{iAn}),100*( upper_ci_temporal{iAn} - mean_values_temporal{iAn}), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5,'Color','k');

end

hold on
plot([-100 100], [-100 100],'--k');

ylim([-10, 70])
xlim([-10, 40])

%% 
%% WITH SEM
colors=jet(13);
colors_transparency=colors;
colors_transparency=colors_transparency;
figure(2)
clf
for iAn=1:length(trials_temporal_mean)
    data_temporal=trials_temporal_mean{iAn};
    data_nasal=trials_nasal_mean{iAn};
    hold on
    %plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    
    %plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', colors_transparency(iAn,:),'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    s1=scatter(100*data_temporal(:),100*data_nasal(:),'o','MarkerFaceColor',colors(iAn,:),'MarkerEdgeColor','k');
    s1.MarkerFaceAlpha = .5;
    

    % SEM
    %errorbar(100*mean_values_temporal{iAn}, 100*(mean_values_nasal{iAn} - sem_values_nasal{iAn}),100*(mean_values_nasal{iAn} + sem_values_nasal{iAn}),100*mean_values_nasal{iAn}, 100*(mean_values_temporal{iAn} - sem_values_temporal{iAn}),100*( mean_values_temporal{iAn} + sem_values_temporal{iAn}), 'o', 'MarkerFaceColor', colors(iAn,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5,'Color','k', 'MarkerSize', 12);
    hold on;

    % Create a plot with horizontal  and vertical CI
    errorbar(100*mean_values_temporal{iAn}, 100*mean_values_nasal{iAn},100*(mean_values_nasal{iAn} - lower_ci_nasal{iAn}),100*( upper_ci_nasal{iAn} - mean_values_nasal{iAn}), 100*(mean_values_temporal{iAn} - lower_ci_temporal{iAn}),100*( upper_ci_temporal{iAn} - mean_values_temporal{iAn}),'o', 'MarkerFaceColor', colors(iAn,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5,'Color','k', 'MarkerSize', 12);

end

hold on
plot([-100 100], [-100 100],'--k');

ylim([-10, 70])
xlim([-10, 40])



%%
% plot
% Create a scatter plot with error bars
figure
x = 1:size(data_example_mean, 1);
errorbar(x, mean_values, mean_values - lower_ci, upper_ci - mean_values, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);

   
%% scatter plot COLOR CODED ANIMAL

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


[pval, h0, stats]=signrank(for_test_nasal_data(:), for_test_temporal_data(:),'Tail','right');
pval

hold on;
plot(mean_all_temporal, mean_all_nasal, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

ylabel('Pupil size change (%) NASAL', 'FontSize',18,'Color','k');
xlabel('Pupil size change (%) TEMPORAL','FontSize',18,'Color','k');
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];
%print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_sessions.pdf']);