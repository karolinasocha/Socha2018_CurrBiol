%% plot pupil response median vs fraction of running
% G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\recalculation_pupil_median\figures\area_dynamics_behavior_condidtions
% pupil_alltrials_alltrials_fractionrunning_average_pupil.pdf

clear all

load('G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour\data_behav_40expt.mat')

newdir=data_behav_40expt.newdir;
stimulus=data_behav_40expt.stimulus;
velocity=data_behav_40expt.velocity;
expt=data_behav_40expt.expt;
expt2=data_behav_40expt.expt2;
clear data_behav_40expt

for iAn=1:size(newdir,1)
clear filepathanalysis
clear behav_name
%clear diameter_data_recalculated_revision
clear diameter_data_recalculated_revision_noartifacts
    
    filepathanalysis=['J:\data\analysis\' newdir{iAn}]; 
    %behav_name = [filepathanalysis, '\diameter_data_recalculated_revision.mat'];
    behav_name = [filepathanalysis, '\diameter_data_recalculated_revision_noartifacts.mat'];

    load(behav_name)
    diameter_data{iAn}=diameter_data_recalculated_revision_noartifacts;
end

for iAn=1:size(newdir,1)
diam_epochs_area{iAn}=diameter_data{iAn}.diameter_area.ddat_interp_area;
diam_epochs_zscore{iAn}=diameter_data{iAn}.diam_zscore.ddat_interp_zscore;
diam_epochs_normmedian{iAn}=diameter_data{iAn}.diam_norm_median.ddat_interp_norm_median;

end
%%
diameter_array=diam_epochs_normmedian;
%% PUPIL INTERPOLATION
clear array_X_data
for iAn=1:size(diameter_array,2);

clear pupil_epoch
clear tmp_template
clear tmp
clear dtime
clear p2time
pupil_epoch=diameter_array{iAn};
tmp=nan(size(pupil_epoch));
tmp_template=nan(315,size(tmp,2),size(tmp,3));

for itrials=1:size(tmp,2)
    for istim=1:size(tmp,3) 
        
    dtime = linspace(0,1,length(pupil_epoch));
    p2time =  linspace(0,1,length(tmp_template));
    tmp_template(:,itrials,istim)= interp1(dtime,pupil_epoch(:,itrials,istim),p2time,'linear','extrap');
    end
end

array_X_data{iAn}=tmp_template;

end
array_X_data_trials=array_X_data %reordered depending on stimulus
%% median change

for iAn=1:40;
av_pupil_median{iAn}=squeeze(nanmean(array_X_data_trials{iAn}(63:223,:,:)));
end

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
    ntrial=size(av_median{iAn},1);
    for iitrial=1:ntrial
    clear tmp_trace
    clear tmp_fraction
    tmp_trace=av_median{iAn}(iitrial,1:end-1);
    tmp_fraction=frac_running_pertrial{iAn}(iitrial);
    tmp_mean_velocity=mean_running_epoch{iAn}(iitrial);

    j=j+1;
    array_fraction_running(j)=tmp_fraction;
    array_av_pupil(j,:)=tmp_trace;
    array_mean_running(j)=tmp_mean_velocity;
    end    
end

%%

  data_set_pupil=cell2mat(av_median');
  data_set_pupil=data_set_pupil(:,1:end-1);
  data_set_running=cell2mat(frac_running_pertrial);
%%
  figure
  clf
    
    ax_heatmap_trials=axes('Position',[0.1 0.2 0.1 0.7]);
    tmp=data_set_pupil;
    [val ord]=sort(data_set_running);
    array_mean_running_sorted=array_mean_running(ord);
    tt = imagesc(tmp(ord,:),[-.5 .5]);
    ylabel('Trials #')
    xlabel('Stimulus (deg)')
    title(sprintf('Pupil norm all sessions'));
    colormap(bluered)
    set(gca,'box','off','tickdir','out','xtick',1:12,'xticklabel',[0:30:330],'XTickLabelRotation',45,...
         'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
   
    colorbar('peer',ax_heatmap_trials,'position',[0.23 0.2 0.02 0.1])
  
    ax_histogram_velocityfraction=axes('Position',[0.30 0.2 0.05 0.7])
    plot(100*val',1:size(val,2),'.r');
    title('Fraction running')     
    xlabel('Fraction running (%)');
    set(gca, 'YDir','reverse');
    hold on
    plot([25 25],[0 size(tmp,1)],'--k');
    plot([75 75],[0 size(tmp,1)],'--k');

    set(gca,'xtick',[0:25:100],'yticklabel',[]);
    axis tight
       set(gca,'box','off','tickdir','out',...
         'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
    %ylabel('#Trials')
%    ax_animalid=axes('Position',[0.450 0.2 0.05 0.7])

%    color_animal_id=jet(13);
%     for indx=1:398
%     c1_an=data_animal_id(ord(indx));
%     plot([0 1],[indx indx],'-','color',color_animal_id(c1_an,:),'linewidth',5);
%     hold on
%     end
%     xlabel('Animal');   
%     set(gca, 'YDir','reverse','xtick',[],'yticklabel',[]);
%     axis tight
%     set(gca,'box','off','tickdir','out',...
%          'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
%     title('Animal')
%     
    
   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure2\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'pupil_alltrials_alltrials_fractionrunning_animalID.pdf']);

   %% plot averages of pupil size 
    tmp_average_trace=data_set_pupil;
    [val ord]=sort(data_set_running);
    tmp_average_trace=tmp_average_trace(ord,:);
    %% plot every 10
  size_bin=50;
  clear av_running
  clear av_trace
  start_epoch=1:50:358
  end_epoch=[start_epoch(2:end)-1 size(tmp_average_trace,1)]
    for nround=1:8
        mean_velocity(nround)=nanmean(array_mean_running_sorted(start_epoch(nround):end_epoch(nround)));
        nanstd_mean_velocity(nround)=nanstd(array_mean_running_sorted(start_epoch(nround):end_epoch(nround))./sqrt(size_bin));
        av_trace(nround,:)=nanmean(tmp_average_trace(start_epoch(nround):end_epoch(nround),:),1);
        av_fraction_running(nround,:)=nanmean(val(start_epoch(nround):end_epoch(nround)));
    end
    
%ax_average_trace=axes('Position',[0.450 0.2 0.05 0.7])
av_mean_running_flip=flip(mean_velocity',1);
nanstd_mean_running_flip=flip(nanstd_mean_velocity',1);
av_fraction_running_flip=flip(av_fraction_running,1);
av_trace_flip=flip(av_trace,1);

height_plot=0.7;
height_subplot=height_plot/nround;
%%
 %ax_average_trace=axes('Position',[0.450 0.2+height_subplot*iiplot 0.05 height_subplot])

 for iiplot=1:8
     ax_average_trace(iiplot)=axes('Position',[0.450 0.2+height_subplot*(iiplot-1) 0.05 height_subplot]);
     hold on
     plot(av_trace_flip(iiplot,:),'k','linewidth',2);
     axis tight
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
    str_fraction_vel=sprintf('fract=%.1f',100*av_fraction_running_flip(iiplot));
    str_mean_vel=sprintf('\nmean=%.2f cm/s',av_mean_running_flip(iiplot));
    str_nanstdf_vel=sprintf('\nstd=%.2f cm/s',nanstd_mean_running_flip(iiplot));
    text(20,-0.5+height_subplot*(iiplot-1),[str_fraction_vel str_mean_vel str_nanstdf_vel]);

 end
title('Av Pupil')     

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure2\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'pupil_heatmaps_pupil_fractionrunning_pupil.pdf']);



