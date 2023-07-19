% Figure 1H scatter plot

%%
clear all
% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')

%% pvals for individual session
clear h0
clear pval
clear ci
clear stats

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

number_sessions=size(new_pupil_data.trials_raw_diam_diff_relative_stimulation,2)

[vals, indx]=sort(new_pupil_data.animal_id);

%% pvals for individual session
nColors=40
colors_scatter=jet(nColors);
figure
for iSess=1:number_sessions
    
    %data_sessions_in_sorted{iSess}=new_pupil_data.trials_raw_diam_diff_relative_stimulation{indx(iSess)};

    clear appendedArray;
    %appendedArray=data_sessions_in_sorted{iSess}
    appendedArray=new_pupil_data.av_raw_diam_diff_relative_stimulation{iSess};
    
    clear data_nasal
    clear data_temporal
    
    data_nasal=appendedArray(:,order_nasal);
    data_temporal=appendedArray(:,order_temporal);
    hold on
    plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    
    data_nasal_all{iSess}=data_nasal;
    data_temporal_all{iSess}=data_temporal;
   
end
hold on
plot([-80 80], [-80 80],'--k');

ylim([-40, 80])
xlim([-40, 40])

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
print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_sessions.pdf']);


%%
