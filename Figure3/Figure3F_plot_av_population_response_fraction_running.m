clear all
analys_path='G:\mousebox\analysis\';

newdir={'140623_KS092_2P_KS\run01_ori_ds_V1',...
    '140623_KS093_2P_KS\run01_ori_ds_V1',...
    '150420_KS133_2P_KS\run01_ori_ds_V1',...
    '150915_KS145_2P_KS\run01_ori_ds_V1_full',...
    '140703_KS099_2P_KS\run01_ori_ds_V1',...
    ...%'141209_KS126_2P_KS\run01_ori_ds_V1_reversed',... %remove for heatmaps  has 9 stims problematic
    '140810_KS103_2P_KS\run01_ori_ds_V1_full',...
    ...%'160105_KS154_2P_KS\run01_orids_V1',... %remove for heatmaps 9 stims
    '151007_KS145_2P_KS\run012_ori_ds_V1_full',...dodane
  ...%'140808_KS092_2P_KS\run02_ori_ds_V1_full',... % remove for heatmaps 9 stims
  '140808_KS093_2P_KS\run02_ori_ds_V1_full',...
  '150420_KS133_2P_KS\run02_ori_ds_V1',...
  '150917_KS145_2P_KS\run02_ori_ds_V1_full',...
  ...%'141215_KS126_2P_KS\run02_ori_ds_V1_reversed',... % remove for heatmaps 9 stims
  '140810_KS103_2P_KS\run02_ori_ds_V1_full',...
  '160219_KS154_2P_KS\run02_ori12_rand_V1',....
  '160223_KS160_2P_KS\run02_ori12_rand_V1',...
  '151007_KS145_2P_KS\run023_ori_ds_V1_full',... % dodane
 '140808_KS092_2P_KS\run03_ori_ds_V1_full', ...
             '140808_KS093_2P_KS\run03_ori_ds_V1_full',  ...
             '150420_KS133_2P_KS\run03_ori_ds_V1',...
             '150915_KS145_2P_KS\run034_ori_ds_V1_full',...
             '140810_KS103_2P_KS\run03_ori_ds_V1_full',...
             '160209_KS154_2P_KS\run03_ori12_rand_V1',...
             '160223_KS160_2P_KS\run03_ori12_rand_V1',...
             '160621_KS166_2P_KS\run03_ori12_V1_awake',...
             '160707_KS164_2P_KS\run03_ori12_V1_awake',...
             '160712_KS167_2P_KS\run03_ori12_V1_awake',...
             ...%'170107_KS174_2P_KS\run03_ori12_V1',... % remove for heatmaps zmienic potem na 170110
             '170110_KS173_2P_KS\run03_ori12_V1_awake',...
             '151007_KS145_2P_KS\run034_ori_ds_V1_full',... % dodane
             '151007_KS145_2P_KS\run03_ori_ds_V1_full',...% dodane
             '170107_KS174_2P_KS\run03_ori12_V1',...% dodane
             '170110_KS174_2P_KS\run03_ori12_V1_awake',...% dodane
             '170108_KS174_2P_KS\run03_ori12_V1_awake',... % dodane
             '170104_KS173_2P_KS\run03_ori12_V1',... % dodane  remove for heatmaps 100 frames
             '170110_KS173_2P_KS\run03_ori12_V1_awake'};% dodane
%%         
for iAn=1:size(newdir,2)
expt{iAn} = frGetExpt(newdir{iAn});
expt2{iAn} =doLoadStimLogs3(expt{iAn});
stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
vel{iAn} = load([expt{iAn}.dirs.analrootpn,'\','velo_final.mat']);
velocity{iAn} = vel{iAn}.velo_final(:,2);
epochs_awake{iAn} = expt2{iAn}.frames.epochs;
stims_awake{iAn} = expt2{iAn}.frames.stims;
     if exist([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat'])
    tcs_awake{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat']);
    tcs_awake{iAn}=tcs_awake{iAn}.tcs_handseg;
    else
         tcs_awake{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_colormap.mat']);
         tcs_awake{iAn}=tcs_awake{iAn}.tcs_colormap;
     end
end  
%%
for iAn=1:length(newdir)
    newdir_animals_id{iAn}=newdir{iAn}(8:12);
end
animal_id_list=unique(newdir_animals_id);
animal_id=nan(length(newdir),1);

for iAn=1:length(animal_id_list)
    indexes=find(strcmp(newdir_animals_id, animal_id_list{iAn}));
    animal_id(indexes)=iAn;
end

%% --------------------------------------

%%
for iAn=1:size(newdir,2)
%tc_anesthesia_zscore{iAn}=zscore(tcss_anesthesia{iAn}.ratio_vis);
tc_awake_ratio{iAn}=lowpass(tcs_awake{iAn}.ratio_vis);
%tc_awake_ratio{iAn}=zscore(lowpass(tcs_awake{iAn}.ratio_vis));

end
%% REORDERING AND INTERPOLATION CALCIUM DATA
for iAn=1:size(newdir,2)
    clear calcium_trace
    clear velocity_trace
    % calculate calcium trace
%calcium_trace=tc_awake_zscore{iAn};
calcium_trace=tc_awake_ratio{iAn};
velocity_trace=velocity{iAn};
% epochs calculate epochs so they match as these for general analysis
[cdat_epochs cdat_stims]=calc_calciumepochs(expt2{iAn},calcium_trace);
[vdat_epochs vdat_stims]=calc_calciumepochs(expt2{iAn},velocity_trace);
% reorder 
clear stimulus
stimulus=[expt2{iAn}.prot.pars.ori]
clear stims_value
stims_value=stimulus(1:end-1);
correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];
%
for i=1:12
  order(i)= find(stims_value==correct_value(i));
  order(13)=13;
end

% cdat_ordered_epochs{iAn}=cdat_stims(:,:,:,order);
% vdat_ordered_epochs{iAn}=vdat_stims(:,:,:,order);
% vdat_ordered_stims{iAn}=vdat_stims(:,:,:,order);


cdat_ordered_epochs{iAn}=cdat_epochs(:,:,:,order);
vdat_ordered_epochs{iAn}=vdat_epochs(:,:,:,order);
vdat_ordered_stims{iAn}=vdat_stims(:,:,:,order);

clear cdat_tmp_interp
cdat_tmp_interp=cdat_ordered_epochs{iAn};
clear interpolated_tmp_epochs
length_epochs=315;
% interpolate
[interpolated_tmp_epochs]=interpolate_calciumepochs(cdat_tmp_interp,length_epochs);
cdat_interpolated_epochs{iAn}=interpolated_tmp_epochs;

clear vdat_tmp_interp
vdat_tmp_interp=vdat_ordered_epochs{iAn};
clear interpolated_tmp_epochs
length_epochs=315;
% interpolate
[v_interpolated_tmp_epochs]=interpolate_calciumepochs(vdat_tmp_interp,length_epochs);
vdat_interpolated_epochs{iAn}=v_interpolated_tmp_epochs;

end

%% FRACTION VELOCITY 

clear frac_running
clear frac_running_pertrial
clear mean_running_epoch
j=0;

for iAn=1:size(newdir,2)
    clear X_velo_tmp
    clear XX_tmp
    XX_tmp=squeeze(vdat_ordered_stims{iAn});
    %XX_tmp_stims=squeeze(vdat_ordered_epochs{iAn});
    sel_iAn=animal_id(iAn);
    for ntrials=1:size(XX_tmp,2)
        j=j+1;
        X_velo_tmp_trial=XX_tmp(:,ntrials,:);
        frac_running_pertrial{iAn}(ntrials)= sum(X_velo_tmp_trial(:)>1,1)/size(X_velo_tmp_trial(:),1);
        mean_running_epoch{iAn}(ntrials)=nanmean(X_velo_tmp_trial(:));

        for iStim=1:size(XX_tmp,3)
            X_velo_tmp=XX_tmp(:,ntrials,iStim);
            %frac_running = sum(X_velo_tmp>1,1)/size(X_velo_tmp,1);
            frac_running{iAn}(ntrials,iStim) = sum(X_velo_tmp(:)>1,1)/size(X_velo_tmp(:),1);
        end
        animal_id_trials(j)=sel_iAn;
    end    
end
%% population plot 
clear cdat_permuted_epochs
clear vdat_permuted_epochs
for iAn=1:size(newdir,2) 
cdat_permuted_epochs{iAn}=permute(cdat_interpolated_epochs{iAn}(:,:,:,1:end-1),[3,1,2,4]);
vdat_permuted_epochs{iAn}=permute(vdat_interpolated_epochs{iAn}(:,:,:,1:end-1),[3,1,2,4]);

end

for iAn=1:size(newdir,2);

    for itrials=1:size(frac_running_pertrial{iAn},2);
    array_tmp{iAn}(:,itrials)=squeeze(nanmean(cdat_permuted_epochs{iAn}(itrials,:,:),2));
    array_vtmp{iAn}(:,itrials)=squeeze(nanmean(vdat_permuted_epochs{iAn}(itrials,:,:),2));

    end
end

frac_running_array=cell2mat(frac_running_pertrial);
calcium_response_array=cell2mat(array_tmp);
velocity_response_array=cell2mat(array_vtmp);
epochlength=size(calcium_response_array,1);
stim_start=2*315/10;
epoch_duration=10*315/10
[val ord]=sort(frac_running_array)
calcium_response_array_ord=calcium_response_array(:,ord);
velocity_response_array_ord=velocity_response_array(:,ord);
animal_id_trials_sorted=animal_id_trials(ord);
%%
clear indx_highlocomotion
clear indx_lowrunning
clear indx_stationary
indx_highlocomotion=find(val>=.250)
%indx_medianrunning=find(val>=0.5 & val<0.75);
indx_lowrunning=find(val>=0.05 & val<0.25);
indx_stationary=find(val>=0 & val<0.05);

%% calculate average + std
% stationary less than 5% running; low 5-25%; high more than 25%
% nanmean(nanmean(velocity_response_array_ord(:,indx_stationary)))
% %0.1439 cm/s
% nanstd(nanmean(velocity_response_array_ord(:,indx_stationary)))
% %0.1903 cm/s
% 
% nanmean(nanmean(velocity_response_array_ord(:,indx_lowrunning)))
% %0.9697 cm/s
% nanstd(nanmean(velocity_response_array_ord(:,indx_lowrunning)))
% %0.6183 cm/s
% 
% 
% nanmean(nanmean(velocity_response_array_ord(:,indx_highlocomotion)))
% %5.3082 cm/s
% nanstd(nanmean(velocity_response_array_ord(:,indx_highlocomotion)))
% %3.7062 cm/s

%%
figure, 
subplot(211)
plot(nanmean(velocity_response_array_ord(:,indx_stationary)',1),'k','linewidth',2);
hold on
plot(nanmean(velocity_response_array_ord(:,indx_lowrunning)',1),'color',[.5 .5 .5],'linewidth',2);
hold on
plot(nanmean(velocity_response_array_ord(:,indx_highlocomotion)',1),'r','linewidth',2);
hold on
xstimulus_start=ones(1,12)*stim_start:epoch_duration:epochlength;
xstimulus_end=ones(1,12)*(stim_start+157.5):epoch_duration:epochlength;

ystimulus_limit=[zeros(1,12);5*ones(1,12)];
ystimulus_limit_top=[1*ones(1,12);1*ones(1,12)];

%plot([xstimulus_start;xstimulus_start],ystimulus_limit,'-b');
%hold on
%plot([xstimulus_end;xstimulus_end],ystimulus_limit,'k');
hold on
plot([xstimulus_start;xstimulus_end],ystimulus_limit_top,'k','linewidth',3);

set(gca,'box','off','tickdir','out','xtick',stim_start:epoch_duration:epochlength,'xticklabel',[0:30:330],'XTickLabelRotation',45,...
    'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
ylabel('Velocity (cm/s)');
 
subplot(212)
plot(nanmean(calcium_response_array_ord(:,indx_stationary)',1),'k','linewidth',2);
hold on
plot(nanmean(calcium_response_array_ord(:,indx_lowrunning)',1),'color',[.5 .5 .5],'linewidth',2);
hold on
plot(nanmean(calcium_response_array_ord(:,indx_highlocomotion)',1),'r','linewidth',2);
hold on
xstimulus_start=ones(1,12)*stim_start:epoch_duration:epochlength;
xstimulus_end=ones(1,12)*(stim_start+157.5):epoch_duration:epochlength;

ystimulus_limit=[zeros(1,12);5*ones(1,12)];
ystimulus_limit_top=[150*ones(1,12);150*ones(1,12)];

%plot([xstimulus_start;xstimulus_start],ystimulus_limit,'-b');
%hold on
%plot([xstimulus_end;xstimulus_end],ystimulus_limit,'k');

hold on
plot([xstimulus_start;xstimulus_end],ystimulus_limit_top,'k','linewidth',3);

set(gca,'box','off','tickdir','out','xtick',stim_start:epoch_duration:epochlength,'xticklabel',[0:30:330],'XTickLabelRotation',45,...
    'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
ylabel('dF/F (%)');

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\dFF_awake_anesthesia_fractionrunning\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'plot_avdFF_velocity_calcium_response_population_awake_allexperiments_topaper.pdf']);

%%

av_population_response_ordered=calcium_response_array_ord;

%% plot average over 50 trials after sorting regarding to fraction of running
colors_line=parula(12)
figure

    ax_histogram_velocityfraction=axes('Position',[0.6 0.2 0.05 0.3])
    plot(100*val',1:size(val,2),'.r');
    title('Fraction running')     
    xlabel('Fraction running (%)');
    set(gca, 'YDir','reverse');
    hold on
    plot([25 25],[0 size(val,2)],'--k');
    plot([75 75],[0 size(val,2)],'--k');
    axis tight
    set(ax_histogram_velocityfraction,'xlim',[0 100],'xtick',[0:25:100],'ytick',[0:50:309]);
    ylim([0 309]);
       set(gca,'box','off','tickdir','out',...
         'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)

for iStim =1:6
    ax_averageresponse_trials(iStim)=axes('Position',[0.32 0.45-(iStim-1)*(0.3/6) 0.2 0.3/6]);
end

start_number=1:50:309;
end_number=50:50:309;
end_number(6)=309;

for iStim =1:6
    
    axes(ax_averageresponse_trials(iStim))
    %plot(nanmean(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1),'k','linewidth',2);
    number_animals=length(unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim))));
    meanData=nanmean(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1);
    %semData=nanstd(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1)./sqrt(size(calcium_response_array_ord(:,start_number(iStim):end_number(iStim)),2));
    semData=nanstd(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1)./sqrt(number_animals);
    
    x = 1:size(meanData,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;
    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    ylim([0 120]);
    hold on
    plot([0 size(calcium_response_array_ord,1)],[50 50],'--k');
    axis off
end

axes(ax_averageresponse_trials(6));
hold on
plot([size(calcium_response_array_ord,1)+20 size(calcium_response_array_ord,1)+20],[0 50],'linewidth',2,'color','k')

axes(ax_averageresponse_trials(1));
xstimulus_start=ones(1,12)*stim_start:epoch_duration:epochlength;
xstimulus_end=ones(1,12)*(stim_start+159):epoch_duration:epochlength;

ystimulus_limit=[zeros(1,12);5*ones(1,12)];
ystimulus_limit_top=[120*ones(1,12);120*ones(1,12)];

%plot([xstimulus_start;xstimulus_start],ystimulus_limit,'-b');
%hold on
%plot([xstimulus_end;xstimulus_end],ystimulus_limit,'k');

hold on
plot([xstimulus_start;xstimulus_end],ystimulus_limit_top,'k','linewidth',3);


set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\dFF_awake_anesthesia_fractionrunning\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'plot_avdFF_velocity_calcium_response_population_awake_allexperiments_fractionrunning_v2_withstimulus.pdf']);

%%
colors_line=parula(12)
figure

    ax_histogram_velocityfraction=axes('Position',[0.6 0.2 0.05 0.3])
    plot(100*val',1:size(val,2),'.r');
    title('Fraction running')     
    xlabel('Fraction running (%)');
    set(gca, 'YDir','reverse');
    hold on
    plot([25 25],[0 size(val,2)],'--k');
    plot([75 75],[0 size(val,2)],'--k');
    axis tight
    set(ax_histogram_velocityfraction,'xlim',[0 100],'xtick',[0:25:100],'ytick',[0:50:309]);
    ylim([0 309]);
       set(gca,'box','off','tickdir','out',...
         'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)
   set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
   print(gcf,'-dpdf',[filepathanalysis, 'Figure3Fright_velocity_fraction_running.pdf']);
  
  %%
  
figure
clf

for iStim =1:6
    ax_averageresponse_trials(iStim)=axes('Position',[0.32 0.45-(iStim-1)*(0.3/6) 0.2 0.3/6]);
end

start_number=1:50:309;
end_number=50:50:309;
end_number(6)=309;

for iStim =1:6
    
    axes(ax_averageresponse_trials(iStim))
    %plot(nanmean(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1),'k','linewidth',2);
    animal_tmp=animal_id_trials_sorted(start_number(iStim):end_number(iStim));
    animal_list_tmp=unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim)));
    clear mean_session_array
    clear mean_animal_sel_calcium
    
    for i=1:length(animal_list_tmp)
        tmp_index=find(animal_tmp==animal_list_tmp(i));
        tmp_calcium=calcium_response_array_ord(:,start_number(iStim):end_number(iStim));
        mean_animal_sel_calcium{i}=nanmean(tmp_calcium(:,tmp_index),2);
    end

    mean_session_array=cell2mat(mean_animal_sel_calcium);
    number_animals=length(unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim))));
    meanData=nanmean(mean_session_array',1);
    %semData=nanstd(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1)./sqrt(size(calcium_response_array_ord(:,start_number(iStim):end_number(iStim)),2));
    semData=nanstd(mean_session_array',1)./sqrt(number_animals);
   
    x = 1:size(meanData,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;
    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    ylim([0 120]);
    hold on
    plot([0 size(calcium_response_array_ord,1)],[50 50],'--k');
    axis off
end

axes(ax_averageresponse_trials(6));
hold on
plot([size(calcium_response_array_ord,1)+20 size(calcium_response_array_ord,1)+20],[0 50],'linewidth',2,'color','k')

axes(ax_averageresponse_trials(1));
xstimulus_start=ones(1,12)*stim_start:epoch_duration:epochlength;
xstimulus_end=ones(1,12)*(stim_start+159):epoch_duration:epochlength;

ystimulus_limit=[zeros(1,12);5*ones(1,12)];
ystimulus_limit_top=[120*ones(1,12);120*ones(1,12)];

%plot([xstimulus_start;xstimulus_start],ystimulus_limit,'-b');
%hold on
%plot([xstimulus_end;xstimulus_end],ystimulus_limit,'k');

hold on
plot([xstimulus_start;xstimulus_end],ystimulus_limit_top,'k','linewidth',3);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
print(gcf,'-dpdf',[filepathanalysis, 'Figure3F_plot_avdFF_calcium_response_population_fractionrunning.pdf']);
print(gcf,'-dpdf',[filepathanalysis, 'Figure3F_plot_avdFF_calcium_response_population_fractionrunning.svg']);


%% plot average over 50 trials after sorting regarding to fraction of running
colors_line=parula(12)
figure

    ax_histogram_velocityfraction=axes('Position',[0.6 0.2 0.05 0.3])
    plot(100*val',1:size(val,2),'.r');
    title('Fraction running')     
    xlabel('Fraction running (%)');
    set(gca, 'YDir','reverse');
    hold on
    plot([25 25],[0 size(val,2)],'--k');
    plot([75 75],[0 size(val,2)],'--k');
    axis tight
    set(ax_histogram_velocityfraction,'xlim',[0 100],'xtick',[0:25:100],'ytick',[0:50:309]);
    ylim([0 309]);
       set(gca,'box','off','tickdir','out',...
         'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3)

for iStim =1:6
    ax_averageresponse_trials(iStim)=axes('Position',[0.32 0.45-(iStim-1)*(0.3/6) 0.2 0.3/6]);
end

start_number=1:50:309;
end_number=50:50:309;
end_number(6)=309;

for iStim =1:6
    
    axes(ax_averageresponse_trials(iStim))
    %plot(nanmean(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1),'k','linewidth',2);
    animal_tmp=animal_id_trials_sorted(start_number(iStim):end_number(iStim));
    animal_list_tmp=unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim)));
    clear mean_session_array
    clear mean_animal_sel_calcium
    
    for i=1:length(animal_list_tmp)
        tmp_index=find(animal_tmp==animal_list_tmp(i));
        tmp_calcium=calcium_response_array_ord(:,start_number(iStim):end_number(iStim));
        mean_animal_sel_calcium{i}=nanmean(tmp_calcium(:,tmp_index),2);
    end

    mean_session_array=cell2mat(mean_animal_sel_calcium);
    number_animals=length(unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim))));
    meanData=nanmean(mean_session_array',1);
    %semData=nanstd(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1)./sqrt(size(calcium_response_array_ord(:,start_number(iStim):end_number(iStim)),2));
    semData=nanstd(mean_session_array',1)./sqrt(number_animals);
   
    x = 1:size(meanData,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;
    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    ylim([0 120]);
    hold on
    plot([0 size(calcium_response_array_ord,1)],[50 50],'--k');
    axis off
end

axes(ax_averageresponse_trials(6));
hold on
plot([size(calcium_response_array_ord,1)+20 size(calcium_response_array_ord,1)+20],[0 50],'linewidth',2,'color','k')

axes(ax_averageresponse_trials(1));
xstimulus_start=ones(1,12)*stim_start:epoch_duration:epochlength;
xstimulus_end=ones(1,12)*(stim_start+159):epoch_duration:epochlength;

ystimulus_limit=[zeros(1,12);5*ones(1,12)];
ystimulus_limit_top=[120*ones(1,12);120*ones(1,12)];

%plot([xstimulus_start;xstimulus_start],ystimulus_limit,'-b');
%hold on
%plot([xstimulus_end;xstimulus_end],ystimulus_limit,'k');

hold on
plot([xstimulus_start;xstimulus_end],ystimulus_limit_top,'k','linewidth',3);

set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
print(gcf,'-dpdf',[filepathanalysis, 'Figure3F_plot_avdFF_velocity_calcium_response_population_awake_allexperiments_fractionrunning_v2_withstimulus.pdf']);

%%

%%
stim_start=63;
stim_end=stim_start+157.5;
    

for iStim =1:6
    
    axes(ax_averageresponse_trials(iStim))
    %plot(nanmean(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1),'k','linewidth',2);
    animal_tmp=animal_id_trials_sorted(start_number(iStim):end_number(iStim));
    animal_list_tmp=unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim)));
    clear mean_session_array
    clear mean_animal_sel_calcium
    
    for i=1:length(animal_list_tmp)
        tmp_index=find(animal_tmp==animal_list_tmp(i));
        tmp_calcium=calcium_response_array_ord(:,start_number(iStim):end_number(iStim));
        mean_animal_sel_calcium{i}=nanmean(tmp_calcium(:,tmp_index),2);
    end

    mean_session_array=cell2mat(mean_animal_sel_calcium);
    number_animals=length(unique(animal_id_trials_sorted(start_number(iStim):end_number(iStim))));
    meanData=nanmean(mean_session_array',1);
    %semData=nanstd(calcium_response_array_ord(:,start_number(iStim):end_number(iStim))',1)./sqrt(size(calcium_response_array_ord(:,start_number(iStim):end_number(iStim)),2));
    semData=nanstd(mean_session_array',1)./sqrt(number_animals);
   
    x = 1:size(meanData,2);

    % Plot the mean data
    plot(x, meanData, 'k', 'LineWidth', 2);
    hold on;
    % Create a shaded region for the SEM
    fill([x, fliplr(x)], [meanData - semData, fliplr(meanData + semData)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    ylim([0 120]);
    hold on
    plot([0 size(calcium_response_array_ord,1)],[50 50],'--k');
    axis off
end


