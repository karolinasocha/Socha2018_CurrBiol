%% create data set with nomalized to median pupil size
%% this is to set analysis on pupil size and assess an effect during stimulation

%% test between animals
%% add directory to functions
addpath(genpath('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\functions'))

%% Main figure 1
clear all
data_behav_40expt_directory='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\data_behav_40expt.mat'
load(data_behav_40expt_directory)
%load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
diameter=data_behav_40expt.diameter;
diam=diameter;
session_id= data_behav_40expt.session_id;
animal_id= data_behav_40expt.animal_id;

%% creates epochs and re-order stimulus so it is consistent for all experiments 

for iAn=1:size(diameter,2)

diam_norm_median=(diam{iAn}-nanmedian(diam{iAn}))./nanmedian(diam{iAn});
diam_substracted_prctile=(diam{iAn}-prctile(diam{iAn},5));
diam_raw=diam{iAn};
diam_zscore=zscore(diam{iAn});

all_epochs=expt2{iAn}.frames.epochs;
stim_epochs=expt2{iAn}.frames.stims;

% raw data
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_raw;
[av_raw_diam_epochs{iAn},trials_raw_diam_epochs{iAn}, av_raw_diam_stims{iAn}, trials_raw_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_raw_diam_epochs_reordered{iAn}, trials_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_raw_diam_epochs{iAn},trials_raw_diam_stims{iAn});
[av_raw_diam_epochs_reordered{iAn}, av_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_raw_diam_epochs{iAn},av_raw_diam_stims{iAn});

% substracted from minimum value
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_substracted_prctile;
[av_prct_diam_epochs{iAn},trials_prct_diam_epochs{iAn}, av_prct_diam_stims{iAn}, trials_prct_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_prct_diam_epochs_reordered{iAn}, trials_prct_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_prct_diam_epochs{iAn},trials_prct_diam_stims{iAn});
[av_prct_diam_epochs_reordered{iAn}, av_prct_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_prct_diam_epochs{iAn},av_prct_diam_stims{iAn});


% norm to median
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_norm_median;
[av_median_diam_epochs{iAn},trials_median_diam_epochs{iAn}, av_median_diam_stims{iAn}, trials_median_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_median_diam_epochs_reordered{iAn}, trials_median_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_median_diam_epochs{iAn},trials_median_diam_stims{iAn});
[av_median_diam_epochs_reordered{iAn}, av_median_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_median_diam_epochs{iAn},av_median_diam_stims{iAn});

end

%% calculate average pupil size during epochs

for iAn=1:length(av_median_diam_epochs_reordered)
    
    [av_median_diam_stims_reordered_onset{iAn}, trials_median_diam_stims_reordered_onset{iAn},...
    av_median_diam_stims_reordered_offset{iAn}, trials_median_diam_stims_reordered_offset{iAn},...
    av_median_diam_stims_reordered_middle{iAn}, trials_median_diam_stims_reordered_middle{iAn},...
    av_median_diam_stims_reordered_entire{iAn}, trials_median_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_median_diam_stims_reordered{iAn}, trials_median_diam_stims_reordered{iAn} );

    [av_prct_diam_stims_reordered_onset{iAn}, trials_prct_diam_stims_reordered_onset{iAn},...
    av_prct_diam_stims_reordered_offset{iAn}, trials_prct_diam_stims_reordered_offset{iAn},...
    av_prct_diam_stims_reordered_middle{iAn}, trials_prct_diam_stims_reordered_middle{iAn},...
    av_prct_diam_stims_reordered_entire{iAn}, trials_prct_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_prct_diam_stims_reordered{iAn}, trials_prct_diam_stims_reordered{iAn} );

    [av_raw_diam_stims_reordered_onset{iAn}, trials_raw_diam_stims_reordered_onset{iAn},...
    av_raw_diam_stims_reordered_offset{iAn}, trials_raw_diam_stims_reordered_offset{iAn},...
    av_raw_diam_stims_reordered_middle{iAn}, trials_raw_diam_stims_reordered_middle{iAn},...
    av_raw_diam_stims_reordered_entire{iAn}, trials_raw_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch(expt2{iAn}, av_raw_diam_stims_reordered{iAn}, trials_raw_diam_stims_reordered{iAn} );

end

%%
for iAn=1:length(av_raw_diam_stims_reordered_entire)
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_raw_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_raw_diam_stims_reordered_offset{iAn};
    trials_raw_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    trials_raw_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    av_raw_diam_diff_relative_stimulation{iAn}=nanmean(trials_raw_diam_diff_relative_stimulation{iAn},1) ;
    
    
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_prct_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_prct_diam_stims_reordered_offset{iAn};
    trials_prct_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    %trials_prct_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    
    clear data_trials_onset
    clear data_trials_offset
    data_trials_onset=trials_median_diam_stims_reordered_onset{iAn};
    data_trials_offset=trials_median_diam_stims_reordered_offset{iAn};
    trials_median_diam_difference_stimulation{iAn}=data_trials_offset - data_trials_onset ;
    %trials_median_diam_diff_relative_stimulation{iAn}=(data_trials_offset - data_trials_onset)./(data_trials_onset) ;
    
end
% %% FIGURE
% for iAn=1:40
% figure(1)
% clf
% subplot(3,1,1)
% plot(squeeze(nanmean(trials_median_diam_difference_stimulation{iAn},1)))
% title('median norm delta')
% 
% subplot(3,1,2)
% plot(squeeze(nanmean(trials_prct_diam_difference_stimulation{iAn},1)))
% title('prct norm delta')
% 
% subplot(3,1,3)
% plot(squeeze(nanmean(trials_raw_diam_difference_stimulation{iAn},1)))
% title('raw delta')
% pause
% end
%%
clear new_pupil_data
new_pupil_data.newdir=data_behav_40expt.newdir;
new_pupil_data.stimulus=data_behav_40expt.stimulus;
new_pupil_data.velocity=data_behav_40expt.velocity;
new_pupil_data.expt=data_behav_40expt.expt;
new_pupil_data.expt2=data_behav_40expt.expt2;
new_pupil_data.diameter=data_behav_40expt.diameter;
new_pupil_data.diam=diameter;
new_pupil_data.sessions_id= data_behav_40expt.session_id;
new_pupil_data.animal_id= data_behav_40expt.animal_id;

new_pupil_data.av_median_diam_stims_reordered_onset=av_median_diam_stims_reordered_onset;
new_pupil_data.trials_median_diam_stims_reordered_onset=trials_median_diam_stims_reordered_onset;
new_pupil_data.av_median_diam_stims_reordered_offset=av_median_diam_stims_reordered_offset;
new_pupil_data.trials_median_diam_stims_reordered_offset=trials_median_diam_stims_reordered_offset;
new_pupil_data.av_median_diam_stims_reordered_middle=av_median_diam_stims_reordered_middle;
new_pupil_data.trials_median_diam_stims_reordered_middle=trials_median_diam_stims_reordered_middle;
new_pupil_data.av_median_diam_stims_reordered_entire=av_median_diam_stims_reordered_entire;
new_pupil_data.trials_median_diam_stims_reordered_entire=trials_median_diam_stims_reordered_entire;

new_pupil_data.trials_median_diam_difference_stimulation=trials_median_diam_difference_stimulation;

new_pupil_data.av_raw_diam_stims_reordered_onset=av_raw_diam_stims_reordered_onset;
new_pupil_data.trials_raw_diam_stims_reordered_onset=trials_raw_diam_stims_reordered_onset;
new_pupil_data.av_raw_diam_stims_reordered_offset=av_raw_diam_stims_reordered_offset;
new_pupil_data.trials_raw_diam_stims_reordered_offset=trials_raw_diam_stims_reordered_offset;
new_pupil_data.av_raw_diam_stims_reordered_middle=av_raw_diam_stims_reordered_middle;
new_pupil_data.trials_raw_diam_stims_reordered_middle=trials_raw_diam_stims_reordered_middle;
new_pupil_data.av_raw_diam_stims_reordered_entire=av_raw_diam_stims_reordered_entire;
new_pupil_data.trials_raw_diam_stims_reordered_entire=trials_raw_diam_stims_reordered_entire;

new_pupil_data.trials_raw_diam_difference_stimulation=trials_raw_diam_difference_stimulation;
new_pupil_data.trials_raw_diam_diff_relative_stimulation=trials_raw_diam_diff_relative_stimulation;
new_pupil_data.av_raw_diam_diff_relative_stimulation=av_raw_diam_diff_relative_stimulation;
    
new_pupil_data.av_prct_diam_stims_reordered_onset=av_prct_diam_stims_reordered_onset;
new_pupil_data.trials_prct_diam_stims_reordered_onset=trials_prct_diam_stims_reordered_onset;
new_pupil_data.av_prct_diam_stims_reordered_offset=av_prct_diam_stims_reordered_offset;
new_pupil_data.trials_prct_diam_stims_reordered_offset=trials_prct_diam_stims_reordered_offset;
new_pupil_data.av_prct_diam_stims_reordered_middle=av_prct_diam_stims_reordered_middle;
new_pupil_data.trials_prct_diam_stims_reordered_middle=trials_prct_diam_stims_reordered_middle;
new_pupil_data.av_prct_diam_stims_reordered_entire=av_prct_diam_stims_reordered_entire;
new_pupil_data.trials_prct_diam_stims_reordered_entire=trials_prct_diam_stims_reordered_entire;

new_pupil_data.trials_prct_diam_difference_stimulation=trials_prct_diam_difference_stimulation;

%save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat','new_pupil_data')

%%

%% FORMAT SO IT FITS FOR STATISTICS

for iAn=1:length(new_pupil_data.trials_raw_diam_diff_relative_stimulation)
    num_trials(iAn)=size(new_pupil_data.trials_raw_diam_diff_relative_stimulation{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_raw_diam_diff_relative_stimulation{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.raw_relative_diam_delta = pupil_nan_array;


%% FORMAT SO IT FITS FOR STATISTICS

for iAn=1:length(new_pupil_data.trials_prct_diam_stims_reordered_offset)
    num_trials(iAn)=size(new_pupil_data.trials_prct_diam_stims_reordered_offset{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_prct_diam_stims_reordered_offset{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.prct_norm_diam_stims_offset = pupil_nan_array;

% ********************************

for iAn=1:length(new_pupil_data.trials_prct_diam_difference_stimulation)
    num_trials(iAn)=size(new_pupil_data.trials_prct_diam_difference_stimulation{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_prct_diam_difference_stimulation{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.prct_norm_diam_stims_delta = pupil_nan_array;

%%

%% SO IT FITS TO STATISTICS
for iAn=1:length(new_pupil_data.trials_raw_diam_stims_reordered_offset);
    num_trials(iAn)=size(new_pupil_data.trials_raw_diam_stims_reordered_offset{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_raw_diam_stims_reordered_offset{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.raw_diam_stims_offset = pupil_nan_array;


% ********************************

% RAW DIFFERENCE

for iAn=1:length(new_pupil_data.trials_raw_diam_difference_stimulation)
    num_trials(iAn)=size(new_pupil_data.trials_raw_diam_difference_stimulation{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_raw_diam_difference_stimulation{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.raw_diam_stims_delta = pupil_nan_array;


%%

%% FOR STATISTICS
for iAn=1:length(new_pupil_data.trials_median_diam_stims_reordered_offset)
    num_trials(iAn)=size(new_pupil_data.trials_median_diam_stims_reordered_offset{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15)

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_median_diam_stims_reordered_offset{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.median_norm_diam_stims_offset = pupil_nan_array;


% ********************************

for iAn=1:length(new_pupil_data.trials_median_diam_difference_stimulation)
    num_trials(iAn)=size(new_pupil_data.trials_median_diam_difference_stimulation{iAn},1);
end

animal_id=new_pupil_data.animal_id;
session_id=new_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=new_pupil_data.trials_median_diam_difference_stimulation{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

new_pupil_data.median_norm_diam_stims_delta = pupil_nan_array;



save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat','new_pupil_data')

% %%
% for iDir=1:12
% average_dir{iDir}=nanmean(new_pupil_data.raw_diam_stims_delta(:,3+iDir));
% end
% 
% [a,b]=signrank(new_pupil_data.raw_diam_stims_delta(:,4),new_pupil_data.raw_diam_stims_delta(:,7), 'Tail', 'right')
% 
% figure, plot(cell2mat(average_dir)*100);
% 
% %%
% 
% for iDir=1:12
% average_dir{iDir}=nanmean(new_pupil_data.raw_relative_diam_delta(:,3+iDir));
% end
% 
% [a,b]=signrank(new_pupil_data.raw_relative_diam_delta(:,4),new_pupil_data.raw_relative_diam_delta(:,7), 'Tail', 'right')
% 
% %%
% %% averaged 40 sessions
% tt=cell2mat(new_pupil_data.av_median_diam_stims_reordered_middle');
% fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
%     'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
% % ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');
% % 
% % axes(ax_heatmap_allmice_plot)
% imagesc(tt)
% axis tight
% axis square
% cmap=colormap(redblue)
% set(gca,'CLim',[-0.2 0.2])
% 
% cb=colorbar;
% cb.Position = [.27,.082,.015,.05] 
% caxis([-0.2 0.3]);
% %%
% 
% %% averaged 40 sessions
% ttrials=cell2mat(new_pupil_data.trials_median_diam_stims_reordered_middle');
% 
% % Find the maximum value in each trial
% max_values = max(ttrials, [], 2);
% 
% % Sort trials based on the maximum value
% [sorted_max_values, sorted_indices] = sort(max_values);
% 
% % Sort the trials matrix based on the sorted indices
% sorted_trials = ttrials(sorted_indices, :);
% 
% %%
% fig(1) = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
%     'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
% % ax_heatmap_allmice_plot = axes('position',[.1,.07,.15,.15],'units','normalized');
% % 
% % axes(ax_heatmap_allmice_plot)
% imagesc(sorted_trials)
% axis tight
% axis square
% cmap=colormap(redblue)
% set(gca,'CLim',[-0.4 0.4])
% 
% cb=colorbar;
% cb.Position = [.02,.082,.015,.05] 
% caxis([-0.2 0.3]);
% 
% %%
% order_nasal=[11,12,1,2,3,4];
% order_temporal=[5,6,7,8,9,10];
% 
% for itrial=1:length(ttrials)
%     [h0(itrial), pval(itrial), ci{itrial}, stats{itrial}]=ttest2(ttrials(itrial,order_nasal), ttrials(itrial,order_temporal),'Tail','right');
%     %[pval(iAn),h0(iAn),stats{iAn}]=signrank(ttrials(itrial,order_nasal), ttrials(itrial,order_temporal), 'Tail','right');
% 
% end
% number_increased_diam_trials=length(find(pval<0.05));
% number_increased_diam_trials