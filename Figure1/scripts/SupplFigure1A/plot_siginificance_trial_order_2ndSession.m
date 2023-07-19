% statistical test for single trial
% plot heatmaps with statistical significance

clear all
% this load processed data for checking statistical tests
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
session_intrest=2 % 1 - first session of the animal with stimulus

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

%print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals_Trial01_2ndsession.pdf']);

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
%print(gcf,'-dpdf',[filepathanalysis, 'heatmap_raw_relative_delta_pupil_1stsession_#trials_2ndSession.pdf']);

%%
