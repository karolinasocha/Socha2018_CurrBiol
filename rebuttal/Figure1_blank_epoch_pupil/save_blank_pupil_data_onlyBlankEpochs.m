%% save data with blank (pupil size)

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
blank_epochs=expt2{iAn}.frames.blanks;

% raw dataames.b
clear eye_diam_entire_trace
eye_diam_entire_trace=diam_raw;
[av_raw_diam_epochs{iAn},trials_raw_diam_epochs{iAn}, av_raw_diam_stims{iAn}, trials_raw_diam_stims{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, stim_epochs);
[trials_raw_diam_epochs_reordered{iAn}, trials_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_raw_diam_epochs{iAn},trials_raw_diam_stims{iAn});
[av_raw_diam_epochs_reordered{iAn}, av_raw_diam_stims_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_raw_diam_epochs{iAn},av_raw_diam_stims{iAn});

[av_raw_diam_epochs{iAn},trials_raw_diam_epochs{iAn}, av_raw_diam_blanks{iAn}, trials_raw_diam_blanks{iAn}]=calculate_average_pupil(eye_diam_entire_trace, all_epochs, blank_epochs);
[trials_raw_diam_epochs_reordered{iAn}, trials_raw_diam_blanks_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, trials_raw_diam_epochs{iAn},trials_raw_diam_blanks{iAn});
[av_raw_diam_epochs_reordered{iAn}, av_raw_diam_blanks_reordered{iAn}]=calculate_same_stimulus_order(expt2{iAn},stimulus{iAn}, av_raw_diam_epochs{iAn},av_raw_diam_blanks{iAn});

end

%%
%% calculate average pupil size during epochs
nsec_starts=2
nsec_ends=4.8
 
for iAn=1:length(av_raw_diam_blanks_reordered)
    
    [av_raw_diam_blanks_reordered_onset{iAn}, trials_raw_diam_blanks_reordered_onset{iAn},...
    av_raw_diam_blanks_reordered_offset{iAn}, trials_raw_diam_blanks_reordered_offset{iAn},...
    av_raw_diam_blanks_reordered_middle{iAn}, trials_raw_diam_blanks_reordered_middle{iAn},...
    av_raw_diam_blanks_reordered_entire{iAn}, trials_raw_diam_blanks_reordered_entire{iAn}]=calculate_average_over_epoch_times(expt2{iAn}, av_raw_diam_blanks_reordered{iAn}, trials_raw_diam_blanks_reordered{iAn}, nsec_starts, nsec_ends);

%     [av_raw_diam_stims_reordered_onset{iAn}, trials_raw_diam_stims_reordered_onset{iAn},...
%     av_raw_diam_stims_reordered_offset{iAn}, trials_raw_diam_stims_reordered_offset{iAn},...
%     av_raw_diam_stims_reordered_middle{iAn}, trials_raw_diam_stims_reordered_middle{iAn},...
%     av_raw_diam_stims_reordered_entire{iAn}, trials_raw_diam_stims_reordered_entire{iAn}]=calculate_average_over_epoch_times(expt2{iAn}, av_raw_diam_stims_reordered{iAn}, trials_raw_diam_stims_reordered{iAn}, nsec_starts, nsec_ends);


end

%%
clear blank_pupil_data
blank_pupil_data.newdir=data_behav_40expt.newdir;
blank_pupil_data.stimulus=data_behav_40expt.stimulus;
blank_pupil_data.velocity=data_behav_40expt.velocity;
blank_pupil_data.expt=data_behav_40expt.expt;
blank_pupil_data.expt2=data_behav_40expt.expt2;
blank_pupil_data.diameter=data_behav_40expt.diameter;
blank_pupil_data.diam=diameter;
blank_pupil_data.sessions_id= data_behav_40expt.session_id;
blank_pupil_data.animal_id= data_behav_40expt.animal_id;

blank_pupil_data.av_raw_diam_blanks_reordered_offset=av_raw_diam_blanks_reordered_offset;
blank_pupil_data.trials_raw_diam_blanks_reordered_offset=trials_raw_diam_blanks_reordered_offset;
blank_pupil_data.trials_raw_diam_blanks_reordered_onset=trials_raw_diam_blanks_reordered_onset;
blank_pupil_data.av_raw_diam_blanks_reordered_onset=av_raw_diam_blanks_reordered_onset;
blank_pupil_data.trials_raw_diam_blanks_reordered_entire=trials_raw_diam_blanks_reordered_entire;
blank_pupil_data.av_raw_diam_blanks_reordered_entire=av_raw_diam_blanks_reordered_entire;
blank_pupil_data.av_raw_diam_blanks_reordered_middle=av_raw_diam_blanks_reordered_middle;
blank_pupil_data.trials_raw_diam_blanks_reordered_middle=trials_raw_diam_blanks_reordered_middle;



%%
% RAW DIFFERENCE

for iAn=1:length(blank_pupil_data.trials_raw_diam_blanks_reordered_offset)
    num_trials(iAn)=size(blank_pupil_data.trials_raw_diam_blanks_reordered_offset{iAn},1);
end

animal_id=blank_pupil_data.animal_id;
session_id=blank_pupil_data.sessions_id;

pupil_nan_array = NaN(sum(num_trials), 15);

t=0

for iAn=1:length(animal_id)
    clear tmp
    clear tmp_control
    
    animal_id_tmp=animal_id(iAn);
    session_id_tmp=session_id(iAn);
    tmp=blank_pupil_data.trials_raw_diam_blanks_reordered_offset{iAn};
    
    for itrial=1:size(tmp,1);
        t=t+1
        pupil_nan_array(t,1)=animal_id_tmp;
        pupil_nan_array(t,2)=session_id_tmp;
        pupil_nan_array(t,3)=itrial;
        pupil_nan_array(t,4:15)=tmp(itrial,:);

    end
end

blank_pupil_data.raw_diam_blanks = pupil_nan_array;
% 
save('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\blank_pupil_data.mat','blank_pupil_data')

new_pupil_data=blank_pupil_data
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
        tmp=new_pupil_data.trials_raw_diam_blanks_reordered_offset{indecies(ii)};
        
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
number_increased_diam_animals


%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS

mean_pupil_delta=cell2mat(pupil_animals_cell');

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

ylabel('Pupil size (mm2)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:30:330],'ylim',[0.2, 1.2]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\Figure1_blank_epoch_pupil\'

print(gcf,'-dpdf',[filepathanalysis, 'Figure1E_average_raw_diam_relative_delta_pupil_SEM_animals_BLANK.pdf']);

%%
%% PLOT MEAN PUPIL CHANGE ACROSS ANIMALS 4 STIMS

mean_pupil_delta=cell2mat(pupil_animals_cell');
mean_pupil_delta=mean_pupil_delta(:,1:3:12); % ranksum ptest  0.3560
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

ylabel('Pupil area (mm2)', 'FontSize',16,'Color','k');
xlabel('Directions (deg)','FontSize',16,'Color','k');
set(gca,'xtick',[1:1:12],'xticklabel',[0:90:330],'ylim',[0.2, 1.2]);
set(gca,'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
% filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure1\scripts\'];

filepathanalysis='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\Figure1_blank_epoch_pupil\'

print(gcf,'-dpdf',[filepathanalysis, 'Figure1E_average_raw_diam_relative_delta_pupil_SEM_animals_BLANK_4stims.pdf']);



