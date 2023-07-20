

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
%% HISTOGRAM PLOT SIGNIFICANT SESSIONS PER ANIMAL
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

%print(gcf,'-dpdf',[filepathanalysis, 'significant_session_per_animal_all.pdf']);

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