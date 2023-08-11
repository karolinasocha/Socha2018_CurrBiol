%% plot differences delta pupil
% Code written 21-07-2023 for Figure 1G
% stationary vs locomotion
%% do calculation of still and running trials
% do reordering of trials so it starts from 0 deg - 330; calculates number
% of running trials depends on stimulus
% calculates pupil diameter for running and non-running trials
% calculates correlation between velocity and pupil diameter
% KS 2018/05/04
% open plot_peak_eye_depends_stim to create trials_eye_stims

%% eye data calculate epochs
% fixed pupil dynamics; before issue because post frames were not taken
% correctly;
clear all

load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')

newdir=new_pupil_data.newdir;
stimulus=new_pupil_data.stimulus;
velocity=new_pupil_data.velocity;
expt=new_pupil_data.expt;
expt2=new_pupil_data.expt2;
diameter=new_pupil_data.diameter;

%% THIS IS ABOUT VELOCITY ONLY
% SELECT TRIALS STILL, LOCOMOTION

% select trials stationary vs locomotion vs non-defined 
% not properly ordered trials

for iAn=1:size(velocity,2);

    vtime=1:length(velocity{iAn});
    %vval=velocity{iAn}(:,2);
    vval=velocity{iAn};
    % parameters to classify the epoch whether still or locomotion
    
    locthresh=1;
    stillthresh=1;
    perthresh=0.95;
    tstart = nan(size(expt2{iAn}.frames.stims));
    tstop = nan(size(expt2{iAn}.frames.stims));
    k=0;
    for i=1:expt2{iAn}.info.nTrials; 
        for j=1:expt2{iAn}.info.nStim;
            k=k+1;
            tstart(i,j)=expt2{iAn}.frames.stims{i,j}(1); %expt2{iAn}.frames.stims{i,j}(1);
            %tstart_vec(k)=expt2{iAn}.frames.stims{i,j}(1)
            tstop(i,j)=expt2{iAn}.frames.stims{i,j}(end);
            %tstop_vec(k)=expt2{iAn}.frames.stims{i,j}(end)
        end; 
    end;

    [loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials(1:length(vval),vval,...
        tstart,tstop,locthresh,stillthresh,perthresh);

end

for iAn=1:size(diameter,2)

velocity_raw=velocity{iAn};

all_epochs=expt2{iAn}.frames.epochs;
stim_epochs=expt2{iAn}.frames.stims;

% here I want to re-order trials for still and locomotion

[trials_still_vel_reordered{iAn}, trials_loc_vels_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, stilltrial{iAn},loctrial{iAn});

end

%%
pupil_data_in=new_pupil_data.trials_raw_diam_diff_relative_stimulation;

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

%% get pupil size per session
pupil_loc_mean=NaN(length(pupil_data_in), 12);
pupil_still_mean=NaN(length(pupil_data_in), 12);
pval_session_loc=nan(length(pupil_data_in),1);
pval_session_still=nan(length(pupil_data_in),1);

for iSess=1:length(pupil_data_in)
    clear pupil_data_tmp
    pupil_data_tmp=pupil_data_in{iSess};
 
    appendedArray_still_tmp= NaN(15, 12);
    appendedArray_loc_tmp=NaN(15, 12);
    
    for iStim=1:expt2{iAn}.info.nStim-1;
    clear tmp_still
    clear tmp_loc
    
    [loc_trial loc_stim_tmp]=find(trials_loc_vels_reordered{iSess}(:,iStim)==1);
    [still_trial still_stim_tmp]=find(trials_still_vel_reordered{iSess}(:,iStim)==1);
    
    appendedArray_still_tmp(1:length(still_trial),iStim)=pupil_data_tmp(still_trial,iStim);
    appendedArray_loc_tmp(1:length(loc_trial),iStim)=pupil_data_tmp(loc_trial,iStim);
    
    toaverage_still=pupil_data_tmp(still_trial,iStim);
    toaverage_loc=pupil_data_tmp(loc_trial,iStim);
    
    pupil_loc_mean(iSess,iStim)=nanmean(toaverage_loc,1);
    pupil_still_mean(iSess,iStim)=nanmean(toaverage_still,1);
    
    end
    
    appendedArray_still{iSess}=appendedArray_still_tmp;
    appendedArray_loc{iSess}=appendedArray_loc_tmp;
    
    clear validIndex
    data_nasal=appendedArray_still{iSess}(:,order_nasal);
    validIndex = ~isnan(data_nasal);
    data_nasal_test=data_nasal(validIndex);
    clear validIndex
    data_temporal=appendedArray_still{iSess}(:,order_temporal);
    validIndex = ~isnan(data_temporal);
    data_temporal_test=data_temporal(validIndex);
    
    %[pval_session(iSess), h0_session(iSess), stats_session{iSess}]=signrank(data_nasal_test(:), data_temporal_test(:),'Tail','right');
    %[h0_session_still(iSess),pval_session_still(iSess),]=ttest2(data_nasal_test(:), data_temporal_test(:),'Tail','right');
     [pval_session_still(iSess), ~, ~]=ranksum(rmmissing(data_nasal_test(:)),rmmissing(data_temporal_test(:)),'Tail','right');

    
    clear validIndex
    data_nasal_loc=appendedArray_loc{iSess}(:,order_nasal);
    validIndex = ~isnan(data_nasal_loc);
    data_nasal_test_loc=data_nasal_loc(validIndex);
    clear validIndex
    data_temporal_loc=appendedArray_loc{iSess}(:,order_temporal);
    validIndex = ~isnan(data_temporal_loc);
    data_temporal_test_loc=data_temporal_loc(validIndex);
    
    %[h0_session_loc(iSess),pval_session_loc(iSess),]=ttest2(data_nasal_test_loc(:), data_temporal_test_loc(:),'Tail','right');
    %[pval_session_loc(iSess), h0_session_loc(iSess), stats_session_loc{iSess}]=signrank(data_nasal_test_loc(:), data_temporal_test_loc(:),'Tail','right');
    try
    [pval_session_loc(iSess), ~, ~]=ranksum(rmmissing(data_nasal_test_loc(:)),rmmissing(data_temporal_test_loc(:)),'Tail','right');
    end
end

%%
for iSess=1:length(appendedArray_loc)
mean_sess_appendedArray_loc(iSess,:)=nanmean(appendedArray_loc{iSess});
mean_sess_appendedArray_still(iSess,:)=nanmean(appendedArray_still{iSess})
end

animal_id_list=unique(new_pupil_data.animal_id);
session_id=new_pupil_data.sessions_id;
animal_id=new_pupil_data.animal_id;

for iAn=1:length(animal_id_list)

    indecies=find(new_pupil_data.animal_id==animal_id_list(iAn));
    
    average_locomotion(iAn,:)=nanmean(mean_sess_appendedArray_loc(indecies,:),1);
    average_stiationary(iAn,:)=nanmean(mean_sess_appendedArray_still(indecies,:),1);
end
%%

mean_all_locomotion=100*nanmean(average_locomotion,1);
sem_all_locomotion=100*nanstd(average_locomotion,1)./sqrt(size(average_locomotion, 1));

mean_all_still=100*nanmean(average_stiationary,1);
sem_all_still=100*nanstd(average_stiationary,1)./sqrt(size(average_stiationary, 1));

x = 1:numel(mean_all_locomotion);

figure(1)
% Plot the mean data
plot(x, mean_all_locomotion, 'r', 'LineWidth', 2);
hold on;
% Create a shaded region for the SEM
fill([x, fliplr(x)], [mean_all_locomotion - sem_all_locomotion, fliplr(mean_all_locomotion + sem_all_locomotion)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, mean_all_still, 'k', 'LineWidth', 2);
hold on;
fill([x, fliplr(x)], [mean_all_still - sem_all_still, fliplr(mean_all_still + sem_all_still)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

axis tight

ylabel('Pupil size change (%)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[-10, 30.1]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

%print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals_stationaryTrials_locomotion.pdf']);

%%
selected_loc_animals=[2,4,11,10]; 
selected_stationary_animals=[2,4,11,10]; 

mean_all_locomotion=100*nanmean(average_locomotion(selected_loc_animals,:),1);
sem_all_locomotion=100*nanstd(average_locomotion(selected_loc_animals,:),1)./sqrt(size(average_locomotion(selected_loc_animals,:), 1));

mean_all_still=100*nanmean(average_stiationary,1);
sem_all_still=100*nanstd(average_stiationary,1)./sqrt(size(average_stiationary, 1));

x = 1:numel(mean_all_locomotion);

figure(1)
% Plot the mean data
plot(x, mean_all_locomotion, 'r', 'LineWidth', 2);
hold on;
% Create a shaded region for the SEM
fill([x, fliplr(x)], [mean_all_locomotion - sem_all_locomotion, fliplr(mean_all_locomotion + sem_all_locomotion)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, mean_all_still, 'k', 'LineWidth', 2);
hold on;
fill([x, fliplr(x)], [mean_all_still - sem_all_still, fliplr(mean_all_still + sem_all_still)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

axis tight

ylabel('Pupil size change (%)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[-10, 30.1]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

%print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals_stationaryAllmice_locomotion4mice.pdf']);

%%
order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

data_locomotion=average_locomotion(selected_loc_animals,:);
data_locomotion_nasal=data_locomotion(:,order_nasal);
data_locomotion_temporal=data_locomotion(:,order_temporal);
%[pval_loc, h0_loc, stats_loc]=signrank(data_locomotion_nasal(:), data_locomotion_temporal(:),'Tail','right');
[pval_loc, h0_loc, stats_loc]=ranksum(rmmissing(data_locomotion_nasal(:)),rmmissing(data_locomotion_temporal(:)),'Tail','right');
% pval_loc= 0.0820 with ranksum
% pval_loc= 0.1236 with signrank
data_stationary=average_stiationary;
data_stationary_nasal=data_stationary(:,order_nasal);
data_stationary_temporal=data_stationary(:,order_temporal);
%[pval_still, h0_still, stats_still]=signrank(data_stationary_nasal(:), data_stationary_temporal(:),'Tail','right');
[pval_still, h0_still, stats_still]=ranksum(rmmissing(data_stationary_nasal(:)),rmmissing(data_stationary_temporal(:)),'Tail','right');
% pval_still= 2.1526e-09 with ranksum
% pval_still= 2.6618e-07 with signrank

%%
animal_id_list=unique(new_pupil_data.animal_id);
y=animal_id_list
% significant_session_index=find(pval_session<0.05);

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

%print(gcf,'-dpdf',[filepathanalysis, 'significant_session_per_animal_stationary.pdf']);
%%
%% average sessions from the same animal

data_nasal=appendedArray(:,order_nasal);
data_temporal=appendedArray(:,order_temporal);

%[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
[pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');

%% AVERAGE EFFECT PER ANIMAL

%% pvals for individual animals
clear h0
clear pval
clear ci
clear stats

order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

animal_id_list=unique(new_pupil_data.animal_id)
for iAn=1:length(animal_id_list)

    indecies=find(new_pupil_data.animal_id==animal_id_list(iAn)); % sessions from same animal
    appendedArray=[];
    
    for ii=1:length(indecies)
        clear tmp
        %tmp=new_pupil_data.trials_median_diam_stims_reordered_offset{indecies(ii)};
        tmp=pupil_still_mean(indecies(ii),:);
        
        appendedArray = vertcat(appendedArray, tmp);
    end
    
    data_nasal=appendedArray(:,order_nasal);
    data_temporal=appendedArray(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    %[pval(iAn), h0(iAn), stats{iAn}]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    [pval(iAn), h0(iAn), stats{iAn}]=ranksum(rmmissing(data_nasal(:)),rmmissing(data_temporal(:)),'Tail','right');
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


%%
%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS

mean_pupil_delta=100*cell2mat(pupil_animals_cell');
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

%print(gcf,'-dpdf',[filepathanalysis, 'average_raw_diam_relative_delta_pupil_SEM_animals_stationaryTrials.pdf']);
%% heatmap

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

%%
figure,

ax1_delta=axes('position',[0.1, 0.2, 0.27,0.16])
shadedErrorBar(1:size(eye_trace_alltrials,2),100*squeeze(nanmean(eye_trace_alltrials,1)),...
    100*nanstd(eye_trace_alltrials,[],1)./sqrt(size(eye_trace_alltrials,1)),'k');
hold on
plot([1 12],[0 0],'--k')
axis tight
hold on
shadedErrorBar(1:size(eye_trace_still,2),100*squeeze(nanmean(eye_trace_still,1)),...
    100*nanstd(eye_trace_still,[],1)./sqrt(size(eye_trace_still,1)),'k');

ylabel('Pupil size change (%)');
set(ax1_delta,'xtick',[1:3:12],'xticklabel',[0:90:330],...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',12,'ticklength',get(ax1_delta,'ticklength')*4);
ylim([-5 10])

box off

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'figure1_deltapupil_all_stationary_40expt.pdf']);

   %% pvals
   order_nasal=[11,12,1,2,3,4];
   order_temporal=[5,6,7,8,9,10];

    data_nasal=eye_trace_still(:,order_nasal);
    data_temporal=eye_trace_still(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    [pval_still, h0_still, stats_still]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    number_increased_diam_animals=length(find(pval_still<0.05));
    number_increased_diam_animals

    data_nasal=eye_trace_alltrials(:,order_nasal);
    data_temporal=eye_trace_alltrials(:,order_temporal);
    
    %[h0(iAn), pval(iAn), ci{iAn}, stats{iAn}]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    %[pval_all, h0_all, stats_all]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    [pval_all, h0_all, stats_all]=ranksum(rmmissing(data_nasal(:)), rmmissing(data_temporal(:)),'Tail','right');
    number_increased_diam_animals=length(find(pval_still<0.05));
    number_increased_diam_animals
    
    %% PLOT PROPORTIONS STATIONARY VS LOCOMOTION
    
%     
% for ii=1:size(stilltrial,2)
% number_still_trials{ii}=sum(stilltrial{ii}(:,1:end-1))
% number_locomotion_trials{ii}=sum(loctrial{ii}(:,1:end-1))
% number_nondefined_trials{ii}=sum(nondefined{ii}(:,1:end-1))
% end
% 
% number_still_trials_array=cell2mat(number_still_trials');
% number_locomotion_trials_array=cell2mat(number_locomotion_trials');
% number_nondefined_trials_array=cell2mat(number_nondefined_trials');