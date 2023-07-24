%% plot pupil response median vs fraction of running
% G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\recalculation_pupil_median\figures\area_dynamics_behavior_condidtions
% pupil_alltrials_alltrials_fractionrunning_average_pupil.pdf

clear all

%load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\data\new_pupil_data.mat')
newdir=new_pupil_data.newdir;
stimulus=new_pupil_data.stimulus;
velocity=new_pupil_data.velocity;
expt=new_pupil_data.expt;
expt2=new_pupil_data.expt2;

diameter_data=new_pupil_data.diameter
diameter_data{10}=diameter_data{10}';

%% it needs delta pupil
trials_raw_diam_diff_relative_stimulation=new_pupil_data.trials_raw_diam_diff_relative_stimulation;
raw_relative_diam_delta=new_pupil_data.raw_relative_diam_delta;
%% VELOCITY
array_data_velo=[];
matrix_data_velo=[];
velocity_data=velocity;

for iAn=1:size(newdir,1);
    clear av
    clear tmp1
    stims{iAn}=expt2{iAn}.frames.stims;  
    [av tmp1]=tcEpochAverage2(velocity_data{iAn},stims{iAn});
    trials_velocity_stims{iAn}=tmp1;
clear tmp_file
clear tmp_file
for iStim=1:size(tmp1,4)
    dtime = linspace(0,1,length(tmp1));
    p2time =  linspace(0,1,159);
tmp_file(:,:,iStim)=interp1(dtime,squeeze(tmp1(:,:,:,iStim)),p2time,'linear');

end
array_data_velo{iAn}=tmp_file;
matrix_data_velo{iAn}=tmp_file(:,:);
end

%% FRACTION RUNNING SINGLE TRIAL 
clear frac_running_pertrial
clear frac_running
clear alltrials_data_eye
clear alltrials_data_eye
alltrials_data_velocity=cell2mat(array_data_velo);
alltrials_data_velocity=alltrials_data_velocity(:,:,1:end-1);
alltrials_data_velocity=alltrials_data_velocity(:,:);
size(alltrials_data_velocity);

for iAn=1:40
    clear X_velo_tmp
    clear XX_tmp
    XX_tmp=array_data_velo{iAn}(:,:,:);

    for ntrials=1:size(XX_tmp,2)
        X_velo_tmp_trial=XX_tmp(:,ntrials,:);
frac_running_pertrial{iAn}(ntrials)= sum(X_velo_tmp_trial(:)>1,1)/size(X_velo_tmp_trial(:),1);
mean_running_epoch{iAn}(ntrials)=nanmean(X_velo_tmp_trial(:));

        for iStim=1:size(XX_tmp,3)
X_velo_tmp=XX_tmp(:,ntrials,iStim);
%frac_running = sum(X_velo_tmp>1,1)/size(X_velo_tmp,1);
frac_running{iAn}(ntrials,iStim) = sum(X_velo_tmp(:)>1,1)/size(X_velo_tmp(:),1);
        end
    end    
end

%%
j=0;
for iAn=1:40    
    ntrial=size(frac_running{iAn},1);
    for iitrial=1:ntrial
    clear tmp_trace
    clear tmp_fraction
    
    tmp_fraction=frac_running_pertrial{iAn}(iitrial);
    tmp_mean_velocity=mean_running_epoch{iAn}(iitrial);

    j=j+1;
    array_fraction_running(j)=tmp_fraction;
    array_mean_running(j)=tmp_mean_velocity;
    end    
end

%%
data_set_pupil=raw_relative_diam_delta;
data_set_running=cell2mat(frac_running_pertrial);

[val ord]=sort(data_set_running);
array_mean_running_sorted=array_mean_running(ord);
data_set_pupil_sorted=data_set_pupil(ord,:);

%%

size_bin=50;
clear av_running
clear av_trace
start_epoch=1:50:358
end_epoch=[start_epoch(2:end)-1 size(data_set_pupil_sorted,1)]

for nround=1:8
    mean_velocity(nround)=nanmean(array_mean_running_sorted(start_epoch(nround):end_epoch(nround)));
    nanstd_mean_velocity(nround)=nanstd(array_mean_running_sorted(start_epoch(nround):end_epoch(nround))./sqrt(size_bin));
    
    mean_pupil(nround,:)=nanmean(data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),4:end),1);
    nanstd_mean_pupil(nround,:)=nanstd(data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),4:end)./sqrt(size_bin));
    clear tmp_pupil_mean
    clear data_nasal
    clear data_temporal
    tmp_pupil_mean=data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),4:end);
    data_nasal=tmp_pupil_mean(:,order_nasal);
    data_temporal=tmp_pupil_mean(:,order_temporal);
    [pval_trials(nround), hval_trials(nround)]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    [hval_trials_tt(nround), pval_trials_tt(nround)]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    
end
%%

    %% plot every 10

height_plot=0.7;
height_subplot=height_plot/nround;
av_trace_flip=flip(mean_pupil,1);
sem_mean_pupil=flip(nanstd_mean_pupil);
%%
 %ax_average_trace=axes('Position',[0.450 0.2+height_subplot*iiplot 0.05 height_subplot])
figure
clf
 for iiplot=1:8
     ax_average_trace(iiplot)=axes('Position',[0.450 0.2+height_subplot*(iiplot-1) 0.05 height_subplot]);
     hold on
     
     semData=sem_mean_pupil(iiplot,:);
     meanData=av_trace_flip(iiplot,:);
     
     x = 1:size(av_trace_flip,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;

    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

     ylim([-0.5 0.5]);
     hold on
     plot([13 13], [-.25 0],'-k','linewidth',3);
     hold on
     plot([1 12],[0 0],'-','color',[0.5 0.5 0.5]);
     axis off
     set(gca, 'xtick',[],'yticklabel',[]);
     set(gca,'box','off','tickdir','out',...
          'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
    set(gca,'Visible','off');

 end
title('Av Pupil')     

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure2\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'Figure2B_pupil_diff_fraction_running.pdf']);

%% average across same animals



order_nasal=[11,12,1,2,3,4];
order_temporal=[5,6,7,8,9,10];

for nround=1:8
    animal_id_sess{nround}=data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),1);
    session_id_sess{nround}=data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),2)';
    animal_names=unique(animal_id_sess{nround});
    clear tmp_pupil_mean
    for iAn=1:length(animal_names)
        indexes=find(animal_id_sess{nround}==animal_names(iAn));
        tmp=data_set_pupil_sorted(start_epoch(nround):end_epoch(nround),4:end);
        tmp_pupil_mean(iAn,:)=nanmean(tmp(indexes,:),1);
    end
    pupil_mean_round{nround}=tmp_pupil_mean;
    data_nasal=tmp_pupil_mean(:,order_nasal);
    data_temporal=tmp_pupil_mean(:,order_temporal);
    [pval(nround), hval(nround)]=signrank(data_nasal(:), data_temporal(:),'Tail','right');
    [hval_tt(nround), pval_tt(nround)]=ttest2(data_nasal(:), data_temporal(:),'Tail','right');
    mean_pupil(nround,:)=nanmean(pupil_mean_round{nround},1);
    nanstd_mean_pupil(nround,:)=nanstd(tmp_pupil_mean./sqrt(length(animal_names)));
    nanimals(nround)=length(animal_names)
    alpha = 0.05;
    alpha_Bonferroni(nround) = alpha / nanimals(nround);
end


av_trace_flip=flip(mean_pupil,1);
sem_mean_pupil=flip(nanstd_mean_pupil);
%%
for nround=1:8

    if pval_tt(nround)<alpha_Bonferroni(nround)
        significance_sign{nround}='*'
    else
        significance_sign{nround}='ns'
    end
end
%significance_sign = {'*'}    {'ns'}    {'*'}    {'*'}    {'*'}    {'*'}    {'ns'}    {'ns'}
%%
ntrials=50
alpha_Bonferroni_trials=0.05/ntrials;
for nround=1:8

    if pval_trials(nround)<alpha_Bonferroni_trials
        significance_sign_trials{nround}='*'
    else
        significance_sign_trials{nround}='ns'
    end
end

%%
 %ax_average_trace=axes('Position',[0.450 0.2+height_subplot*iiplot 0.05 height_subplot])
figure
clf
 for iiplot=1:8
     ax_average_trace(iiplot)=axes('Position',[0.450 0.2+height_subplot*(iiplot-1) 0.05 height_subplot]);
     hold on
     
     semData=sem_mean_pupil(iiplot,:);
     meanData=av_trace_flip(iiplot,:);
     
     x = 1:size(av_trace_flip,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;

    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

     ylim([-0.5 0.5]);
     hold on
     plot([13 13], [-.25 0],'-k','linewidth',3);
     hold on
     plot([1 12],[0 0],'-','color',[0.5 0.5 0.5]);
     axis off
     set(gca, 'xtick',[],'yticklabel',[]);
     set(gca,'box','off','tickdir','out',...
          'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
    set(gca,'Visible','off');
    %str_fraction_vel=sprintf('fract=%.1f',100*av_fraction_running_flip(iiplot));
    %str_mean_vel=sprintf('\nmean=%.2f cm/s',av_mean_running_flip(iiplot));
    %str_nanstdf_vel=sprintf('\nstd=%.2f cm/s',nanstd_mean_running_flip(iiplot));
    %text(20,-0.5+height_subplot*(iiplot-1),[str_fraction_vel str_mean_vel str_nanstdf_vel]);

 end
title('Av Pupil')     

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure2\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'Figure2B_pupil_diff_fraction_running_SEMperAnimal.pdf']);


