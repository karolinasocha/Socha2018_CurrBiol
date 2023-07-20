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
    pupil_animals_cell{iAn}=squeeze(nanmean(appendedArray,1))
    
end

number_increased_diam_animals=length(find(pval<0.05));
number_increased_diam_animals

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

%% across 13 mice
tt=cell2mat(pupil_animals_cell');

[val, idx]=sort(new_pupil_data.animal_id); % sort animals

fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);


imagesc(100*tt)
axis tight
axis square
cmap=colormap(redblue)
set(gca,'CLim',[-0.4 0.4])
ylabel('#Animals', 'FontSize',16)
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
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_raw_relative_delta_pupil_animals.pdf']);
%%
%% scatter plot

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
%    plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
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
%%
order_nasal
order_temporal

colors_scatter=[0.8941    0.1010    0.1098
    0.2157    0.4941    0.7216
    0.3010    0.6863    0.2902
    0.5961    0.3059    0.6392
    1.0000    0.4980    0.0000
    1.0000    1.0000    0.2000];

figure

for iCond=1:12
    pupil_stim{iCond}=av_delta_animals(:,iCond)
end

figure
clf

for iStim=1:6
    clear data_temporal
    clear data_nasal
    
    data_nasal=pupil_stim{order_nasal(iStim)};
    data_temporal=pupil_stim{order_temporal(iStim)};
    hold on
%   plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', [0.7,0.7,0.7],'MarkerEdgeColor', [0.5, 0.5, 0.5]); %colors_scatter(iSess,:)
    plot(100*data_temporal(:),100*data_nasal(:),'o','MarkerSize', 6,'MarkerFaceColor', colors_scatter(iStim,:),'MarkerEdgeColor', [0.5, 0.5, 0.5],); %colors_scatter(iSess,:)

    data_nasal_all{iAn}=nanmean(data_nasal);
    data_temporal_all{iAn}=nanmean(data_temporal);
end

ylim([-10, 70])
xlim([-10, 40])

hold on;
plot(mean_all_temporal, mean_all_nasal, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hLegend=legend
legendLabels = {'300-120', '330-150','0-180','30-210', '60-240','90-270','mean'};
set(hLegend, 'String', legendLabels);

hold on
plot([-100 100], [-100 100],'--k');

ylabel('Pupil size change (%) NASAL', 'FontSize',18,'Color','k');
xlabel('Pupil size change (%) TEMPORAL','FontSize',18,'Color','k');
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];
print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_scatter_13mice.pdf']);
%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS

mean_pupil_delta=100*cell2mat(pupil_animals_cell');

% Calculate the mean and standard deviation across columns
meanData = nanmean(mean_pupil_delta);
stdData = nanstd(mean_pupil_delta);

% Calculate the Standard Error of the Mean (SEM)
semData = stdData ./ sqrt(size(mean_pupil_delta, 1));

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
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals.pdf']);

% calculate mean + SEM
% meanData =[16.1359   12.2291    8.1056    4.2092   -0.3538    1.4475    0.2671   -1.1895   -0.1307    3.9853    4.9183   23.9914]
%semData=[4.1998    2.1429    2.2035    1.7687    1.6019    3.4430    1.3093    1.1727    1.0122    1.8115    1.7518    5.5431 ]

% 0deg: 16.1359 +- 4.1998 (across 13 animals)
% 90deg: 4.2092 +- 1.7687 (across 13 animals)
% 180deg: 0.2671 +- 1.3093 (across 13 animals)
% 270deg: 3.9853 +- 1.8115 (across 13 animals)

%%

%% pvals for individual session
clear h0
clear pval
clear ci
clear stats

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

number_sessions=size(new_pupil_data.trials_raw_diam_diff_relative_stimulation,2)

[vals, indx]=sort(new_pupil_data.animal_id)

for iSess=1:number_sessions
    
    data_sessions_in_sorted{iSess}=new_pupil_data.trials_raw_diam_diff_relative_stimulation{indx(iSess)};

    clear appendedArray;
    appendedArray=data_sessions_in_sorted{iSess}
    
    clear data_nasal
    clear data_temporal
    
    data_nasal=appendedArray(:,order_nasal);
    data_temporal=appendedArray(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    [pval_sess(iSess), h0_sess(iSess), stats_sess{iSess}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    
    pupil_sessions_cell{iSess}=squeeze(nanmean(appendedArray,1))
    
end

number_increased_diam_sessions=length(find(pval_sess<0.05));
number_increased_diam_sessions

%
for iAn=1:length(pval_sess)
    
    if pval_sess(iAn)>0.05
        level_sign_sess{iAn}='ns';
    elseif pval_sess(iAn)<=0.05 &  pval_sess(iAn)>0.01
        level_sign_sess{iAn}='*';
    elseif pval_sess(iAn)<=0.01 &  pval_sess(iAn)>0.001    
        level_sign_sess{iAn}='**';
    elseif pval_sess(iAn)<=0.001 &  pval_sess(iAn)>0.0001    
        level_sign_sess{iAn}='***';
    elseif pval_sess(iAn)<=0.0001    
        level_sign_sess{iAn}='****';
    end

end

%%
%%
y=animal_id_list
significant_session_index=find(pval_sess<0.05);


[x1, y1] =histc(new_pupil_data.animal_id,animal_id_list);
[x2, y2] =histc(new_pupil_data.animal_id(significant_session_index),animal_id_list);

figure, 
bar(animal_id_list,x1, 'FaceColor', [0.5 0.5 0.5]);
hold on
bar(animal_id_list,x2, 'FaceColor', [0.5 0.2 0.8]);

ylabel('#sessions', 'Fontsize', 16);
xlabel('mouse id','Fontsize', 16);

set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

print(gcf,'-dpdf',[filepathanalysis, 'significant_session_per_animal_all.pdf']);

%% PLOT HEAT MAP averaged 40 sessions
%tt=cell2mat(new_pupil_data.av_median_diam_stims_reordered_offset');
[val, idx]=sort(new_pupil_data.animal_id); % sort animals
y_draw_line=[0, find(diff(val)==1)];
x_draw_line=[0.5, 12.5];

tt=cell2mat(pupil_sessions_cell');

fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);

% axes(ax_heatmap_allmice_plot)
imagesc(100*tt(idx,:));
hold on
for ii=1:length(y_draw_line)
plot(x_draw_line,[y_draw_line(ii), y_draw_line(ii)]+0.5,  '-k', 'LineWidth', 1);
hold on
end

ylabels_animal_id={'#1','#2','#3','#4','#5','#6','#7','#8','#9','#10','#11','#12','#13'};
yticks_position=y_draw_line+0.5;
hold on

axis tight
axis square
cmap=colormap(redblue);
set(gca,'CLim',[-20 20]);
ylabel('#Sessions', 'FontSize',16);
xlabel('Directions (deg)','FontSize',16);
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
set(gca,'ytick',yticks_position,'yticklabels',ylabels_animal_id);

% Set tick locations and labels on the right side of the y-axis
% yticksRight = [0:1:40];  % Replace with your own tick locations
% ytickLabelsRight = level_sign_sess;  % Replace with your own tick labels
% set(gca, 'YTick', yticksRight, 'YTickLabel', ytickLabelsRight, 'YAxisLocation', 'right');

cb=colorbar;
cb.Position = [.1,.082,.015,.05];
caxis([-20 20]);

box off

set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);
box off
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

print(gcf,'-dpdf',[filepathanalysis, 'heatmap_raw_diam_relative_delta_pupil_sessions.pdf']);

% data_nasal_sess=tt(:,order_nasal);
% data_temporal_sess=tt(:,order_temporal);
% 
% %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
% [pval_sess_all, h0_sess_all, stats_sess_all]=signrank(data_nasal(:), data_temporal(:),'Tail','right');

    
%% calculate significance per session (take all trials from the single session)
% calculate significance per animal
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];
clear pval

for itrial=1:length(sorted_trials)
    [h0(itrial), pval(itrial), ci{itrial}, stats{itrial}]=ttest2(sorted_trials(itrial,order_nasal), sorted_trials(itrial,order_temporal),'Tail','right');
    %[pval(iAn),h0(iAn),stats{iAn}]=signrank(ttrials(itrial,order_nasal), ttrials(itrial,order_temporal), 'Tail','right');
end

number_increased_diam_sessions=length(find(pval<0.05));
number_increased_diam_sessions
%%
kk=new_pupil_data.trials_median_diam_stims_reordered_offset
clear pval
for isession=1:length(kk)
    clear data_inputs
    clear data_temporal
    clear data_nasal
    data_inputs=kk{isession}
    data_nasal=data_inputs(:,order_nasal)
    data_temporal=data_inputs(:,order_temporal)
    [h0(isession), pval(isession), ci{isession}, stats{isession}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
end
number_increased_diam_sessions=length(find(pval<0.05));
number_increased_diam_sessions

%pval(sorted_indices)
%% averaged 40 sessions
ttrials=cell2mat(new_pupil_data.trials_median_diam_stims_reordered_middle');
% Find the maximum value in each trial
max_values = max(ttrials, [], 2);
% Sort trials based on the maximum value
[sorted_max_values, sorted_indices] = sort(max_values);
% Sort the trials matrix based on the sorted indices
sorted_trials = ttrials(sorted_indices, :);
%
fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
% ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');
% 
% axes(ax_heatmap_allmice_plot)
imagesc(sorted_trials)
axis tight
axis square
cmap=colormap(redblue)
set(gca,'CLim',[-0.4 0.4])
ylabel('# trials', 'FontSize',16)
xlabel('directions','FontSize',16)
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330])
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
cb=colorbar;
cb.Position = [.1,.082,.015,.05] 
caxis([-0.4 0.41]);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\']; 
print(gcf,'-dpdf',[filepathanalysis, 'heatmap_median_pupil_trials.pdf']);


%%
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

for itrial=1:length(ttrials)
    [h0(itrial), pval(itrial), ci{itrial}, stats{itrial}]=ttest2(ttrials(itrial,order_nasal), ttrials(itrial,order_temporal),'Tail','right');
    %[pval(iAn),h0(iAn),stats{iAn}]=signrank(ttrials(itrial,order_nasal), ttrials(itrial,order_temporal), 'Tail','right');

end
number_increased_diam_trials=length(find(pval<0.05));
number_increased_diam_trials