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
    clear indecies
    indecies=find(new_pupil_data.animal_id==animal_id_list(iAn))
    appendedArray=[];
    
    % this takes all trials from multiple sessions
    for ii=1:length(indecies)
        clear tmp
        %tmp=new_pupil_data.trials_median_diam_stims_reordered_offset{indecies(ii)};
        tmp=new_pupil_data.trials_raw_diam_diff_relative_stimulation{indecies(ii)}; % gets all trials
        
        appendedArray = vertcat(appendedArray, tmp);
    end
    
    % this create array with all trials to be compared
    % depending whether stimuli was nasal or temporal directions
    % 
    data_nasal=appendedArray(:,order_nasal);
    data_temporal=appendedArray(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    [pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    pupil_animals_cell{iAn}=squeeze(nanmean(appendedArray,1))
    
    [pval_signrank_animal(iAn), ~, stats_signrank_animal{iAn}]=signrank(100*data_nasal(:), 100*data_temporal(:),'Tail','right');
    
    [~,pval_kstest_animal(iAn), stats_kstest_animal{iAn}]=kstest2(100*rmmissing(data_nasal(:)), 100*rmmissing(data_temporal(:)),'Tail','larger');
    
    [pval_ranksum_animal(iAn), ~, stats_ranksum_animal] = ranksum(100*rmmissing(data_nasal(:)),100*rmmissing(data_temporal(:)),'Tail','right');

    [~, pval_ttest2_animal(iAn)]=ttest2(100*rmmissing(data_nasal(:)),100*rmmissing(data_temporal(:)),'Tail','right');

end

number_increased_diam_animals=length(find(pval<0.05));
number_increased_diam_animals
[pval_signrank_animal; pval_ttest2_animal; pval_ranksum_animal; pval_kstest_animal];

%%
%% scatter plot

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

av_delta_animals=cell2mat(pupil_animals_cell')

for iAn=1:length(pupil_animals_cell)
    clear data_temporal
    clear data_nasal
    
    data_nasal=pupil_animals_cell{iAn}(order_nasal)
    data_temporal=pupil_animals_cell{iAn}(order_temporal)
    hold on
    data_nasal_all{iAn}=nanmean(data_nasal);
    data_temporal_all{iAn}=nanmean(data_temporal);
end

for_test_nasal_data=cell2mat(data_nasal_all');
for_test_temporal_data=cell2mat(data_temporal_all');

% for plots calculate average point 
mean_all_nasal=100*nanmean(for_test_nasal_data(:),1);
sem_all_nasal=100*nanstd(for_test_nasal_data(:),1)./sqrt(size(for_test_nasal_data, 1));

mean_all_temporal=100*nanmean(for_test_temporal_data(:),1);
sem_all_temporal=100*nanstd(for_test_temporal_data(:),1)./sqrt(size(for_test_temporal_data, 1));


%% pvals for individual session
clear h0
clear pval
clear ci
clear stats

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

number_sessions=size(new_pupil_data.trials_raw_diam_diff_relative_stimulation,2)

% sort based on animal ID
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
    [pval_signrank_sess(iSess), ~, stats_signrank_sess{iSess}]=signrank(100*data_nasal(:), 100*data_temporal(:),'Tail','right');
    
    [~,pval_kstest_sess(iSess), stats_kstest_sess{iSess}]=kstest2(100*rmmissing(data_nasal(:)), 100*rmmissing(data_temporal(:)),'Tail','larger');
    
    [pval_ranksum_sess(iSess), ~, stats_ranksum_sess{iSess}] = ranksum(100*rmmissing(data_nasal(:)),100*rmmissing(data_temporal(:)),'Tail','right');

    [~, pval_ttest2_sess(iSess)]=ttest2(100*rmmissing(data_nasal(:)),100*rmmissing(data_temporal(:)),'Tail','right');

    pupil_sessions_cell_sorted{iSess}=squeeze(nanmean(appendedArray,1));
    
end

clear pval_sess
pval_sess=pval_ranksum_sess;

number_increased_diam_sessions=length(find(pval_sess<0.05));
number_increased_diam_sessions

[pval_signrank_sess; pval_ttest2_sess; pval_ranksum_sess; pval_kstest_sess];

%
for iSess=1:number_sessions
    if sum(pval_sess(iSess)>=0.05)==1
    pval_sess_sign{iSess}='ns'
    elseif sum(pval_sess(iSess)<0.05)==1
    pval_sess_sign{iSess}='*'
    end
end


%%
%% PLOT HEAT MAP averaged 40 sessions
%tt=cell2mat(new_pupil_data.av_median_diam_stims_reordered_offset');
[val, idx]=sort(new_pupil_data.animal_id); % sort animals
y_draw_line=[0, find(diff(val)==1)];
x_draw_line=[0.5, 12.5];

tt=cell2mat(pupil_sessions_cell_sorted');
pval_sess_sorted=pval_sess;

clear pval_sess_sign
for iSess=1:number_sessions

    if sum(pval_sess_sorted(iSess)>=0.05)==1
    pval_sess_sign{iSess}='ns'
    elseif sum(pval_sess_sorted(iSess)<0.05)==1
    pval_sess_sign{iSess}='*'
    end
end


fig(2) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
clf
yyaxis right
% axes(ax_heatmap_allmice_plot)
imagesc(100*tt(idx,:));
cmap=colormap(redblue);
set(gca,'CLim',[-20 20]);

hold on
for ii=1:length(y_draw_line)
plot(x_draw_line,[y_draw_line(ii), y_draw_line(ii)]+0.5,  '-k', 'LineWidth', 1);
hold on
end

hold on

% Add labels on the right side of the y-axis
y_tick_vals = 1.0:1:40.5; % Y-axis tick positions
y_tick_labels = pval_sess_sign; % Initialize an empty cell array for y-axis tick labels

set(gca, 'YTick', y_tick_vals, 'YTickLabel', y_tick_labels);

ylabels_animal_id={'#1','#2','#3','#4','#5','#6','#7','#8','#9','#10','#11','#12','#13'};
yticks_position=y_draw_line+0.5;
hold on

yyaxis left
hold on
for ii=1:length(y_draw_line)
plot(x_draw_line,[y_draw_line(ii), y_draw_line(ii)]+0.5,  '-k', 'LineWidth', 1);
hold on
end

hold on
set(gca,'ytick',yticks_position,'yticklabels',ylabels_animal_id);
set(gca, 'Ydir', 'reverse')

axis tight
axis square
ylabel('#Sessions per animal', 'FontSize',16);
xlabel('Directions (deg)','FontSize',16);
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);

cb=colorbar;
cb.Position = [.1,.082,.015,.05];
caxis([-20 20]);

box off

set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);
box off
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\'];
print(gcf,'-dpdf',[filepathanalysis, 'newFigure1D_heatmap_sessions_animal_pval.pdf']);

%%
% 
% %% PLOT HEAT MAP averaged 40 sessions
% %tt=cell2mat(new_pupil_data.av_median_diam_stims_reordered_offset');
% [val, idx]=sort(new_pupil_data.animal_id); % sort animals
% y_draw_line=[0, find(diff(val)==1)];
% x_draw_line=[0.5, 12.5];
% 
% tt=cell2mat(pupil_sessions_cell');
% pval_sess_sorted=pval_sess(idx)
% 
% fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
%     'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
% 
% % axes(ax_heatmap_allmice_plot)
% imagesc(100*tt(idx,:));
% 
% hold on
% for ii=1:length(y_draw_line)
% plot(x_draw_line,[y_draw_line(ii), y_draw_line(ii)]+0.5,  '-k', 'LineWidth', 1);
% hold on
% end
% 
% for iSess=1:number_sessions
%     clear pval_sess_sign
%     if sum(pval_sess(iSess)>=0.05)==1
%     pval_sess_sign='ns'
%     elseif sum(pval_sess(iSess)<0.05)==1
%     pval_sess_sign='*'
%     end
% end
% 
% ylabels_animal_id={'#1','#2','#3','#4','#5','#6','#7','#8','#9','#10','#11','#12','#13'};
% yticks_position=y_draw_line+0.5;
% hold on
% 
% axis tight
% axis square
% cmap=colormap(redblue);
% set(gca,'CLim',[-20 20]);
% ylabel('#Sessions', 'FontSize',16);
% xlabel('Directions (deg)','FontSize',16);
% set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330]);
% set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
% set(gca,'ytick',yticks_position,'yticklabels',ylabels_animal_id);
% 
% % Set tick locations and labels on the right side of the y-axis
% % yticksRight = [0:1:40];  % Replace with your own tick locations
% % ytickLabelsRight = level_sign_sess;  % Replace with your own tick labels
% % set(gca, 'YTick', yticksRight, 'YTickLabel', ytickLabelsRight, 'YAxisLocation', 'right');
% 
% cb=colorbar;
% cb.Position = [.1,.082,.015,.05];
% caxis([-20 20]);
% 
% box off
% 
% set(cb,'tickdir','out','fontsize',8,'ticklength',get(cb,'ticklength')*4);
% box off
% set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
% filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];
% 
% %print(gcf,'-dpdf',[filepathanalysis, 'heatmap_raw_diam_relative_delta_pupil_sessions.pdf']);
% 