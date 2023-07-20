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
