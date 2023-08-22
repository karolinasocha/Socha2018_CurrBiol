%% load data
clear all
filepathanalysis='G:\mousebox\code\mouselab\users\karolina\figure_quantif_behaviour_eye_velocity_selectedexpt';
load([filepathanalysis, '\diameter_data_27expt_v2.mat']);
load([filepathanalysis, '\velocity_27expt_v2.mat']);
load([filepathanalysis, '\expt_27expt_v2.mat']);
load([filepathanalysis, '\expt2_27expt_v2.mat']);
load([filepathanalysis, '\newdir_27expt_v2.mat']);
load([filepathanalysis, '\diameter_norm_data_27expt_v2.mat']);

for iAn=1:23
    
    diameter_data{iAn}=diameter_data_27expt_v2{iAn};
    diameter_norm_data{iAn}=diameter_norm_data_27expt_v2{iAn};
    expt{iAn}=expt_27expt_v2{iAn};
    expt2{iAn}=expt2_27expt_v2{iAn};
    velocity_data{iAn}=velocity_27expt_v2{iAn};
    newdir{iAn}=newdir_27expt_v2{iAn};
    nStim{iAn}=expt2{iAn}.info.nStim;
     if exist([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat'])
    tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat']);
    tcs{iAn}=tcs{iAn}.tcs_handseg;
    else
         tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_colormap.mat']);
         tcs{iAn}=tcs{iAn}.tcs_colormap;
     end
         
end
%%
velocity_data1=velocity_27expt_v2;
for iAn=1:23
    clear an_stat
        an_stat=load([expt{iAn}.dirs.analrootpn,'\stats_version_mean_values.mat']);
        stats_boutons{iAn}=[an_stat.stats_version_mean_values.stats_boutons];
        stats_boutons_normalized{iAn}=[an_stat.stats_version_mean_values.stats_boutons_normalized];
        tc_boutons{iAn}=tcs{iAn}.ratio_vis;
        tc_boutons_zscored{iAn}=zscore(tcs{iAn}.ratio_vis);
        stims{iAn}=expt2{iAn}.frames.stims;
        epochs{iAn}=expt2{iAn}.frames.epochs;
        [tc_response_epochs{iAn} tc_trial_response_epochs{iAn}]=tcEpochAverage2(tc_boutons{iAn},epochs{iAn});
        [tc_response_stims{iAn} tc_trial_response_stims{iAn}]=tcEpochAverage2(tc_boutons{iAn},stims{iAn});
        
        [tc_pupil_stims{iAn} tc_trial_pupil_stims{iAn}]=tcEpochAverage2(diameter_norm_data{iAn},stims{iAn});
        [tc_pupil_epochs{iAn} tc_trial_pupil_epochs{iAn}]=tcEpochAverage2(diameter_norm_data{iAn},epochs{iAn});
        
        [tc_velocity_stims{iAn} tc_trial_velocity_stims{iAn}]=tcEpochAverage2(velocity_data1{iAn},stims{iAn});
        [tc_velocity_epochs{iAn} tc_trial_velocity_epochs{iAn}]=tcEpochAverage2(velocity_data1{iAn},epochs{iAn});

end
%% example experiment check
velocity=velocity_data1;
for iAn=1:23;
diameter_example=diameter_norm_data{iAn};
velocity_example=velocity{iAn};
tc_boutons_zscored_example=lowpass(tc_boutons_zscored{iAn});
tc_boutons_ratio_example=lowpass(tc_boutons{iAn});
tc_population_zscored_example=nanmean(tc_boutons_zscored_example,2);

tc_trial_response_stims_example=tc_trial_response_stims{iAn};
tc_trial_pupil_stims_example=tc_trial_pupil_stims{iAn};
tc_trial_velocity_stims_example=tc_trial_velocity_stims{iAn};

tc_trialav_velocity=nanmean(tc_trial_velocity_stims_example(:,:,:,1:end-1),1);
tc_trialav_pupil=nanmean(tc_trial_pupil_stims_example(:,:,:,1:end-1),1);
tc_trialav_boutons=nanmean(tc_trial_response_stims_example(:,:,:,1:end-1),1);
tc_trialav_population=nanmean(nanmean(tc_trial_response_stims{iAn}(:,:,:,1:end-1),1),2);

clear corr_population_pupil
clear corr_velocity_pupil
clear corr_population_velocity
corr_population_pupil=corrcoef(tc_trialav_population(:), tc_trialav_pupil(:),'rows','complete')
corr_velocity_pupil=corrcoef(tc_trialav_velocity(:), tc_trialav_pupil(:),'rows','complete')
corr_population_velocity=corrcoef(tc_trialav_velocity(:), tc_trialav_population(:),'rows','complete')

clear corr_bouton_pupil
clear corr_bouton_population
clear corr_bouton_velocity

for ib=1:size(tc_trialav_boutons,2)
    clear tmp_r
    clear tmp_r_ex
    tmp_r=squeeze(tc_trialav_boutons(:,ib,:,:));   
    tmp_r_ex=tmp_r';
    corr_bouton_pupil(ib,:,:)=corrcoef(tmp_r(:),tc_trialav_pupil(:),'rows','complete');
    corr_bouton_population(ib,:,:)=corrcoef(tmp_r(:),tc_trialav_population(:),'rows','complete');
    corr_bouton_velocity(ib,:,:)=corrcoef(tmp_r(:),tc_trialav_velocity(:),'rows','complete');
    model_pupil=(nanmean(squeeze(tc_trialav_pupil),1))';
    model_population=(nanmean(squeeze(tc_trialav_population),1))';
    model_speed=(nanmean(squeeze(tc_trialav_velocity),1))';
    
    [ev_bouton_pupil(ib),pv_bouton_pupil(ib)] = calcExpVar_KS(tmp_r_ex,model_pupil);
    [ev_bouton_speed(ib),pv_bouton_speed(ib)] = calcExpVar_KS(tmp_r_ex,model_speed);
    [ev_bouton_population(ib),pv_bouton_population(ib)] = calcExpVar_KS(tmp_r_ex,model_population);
end
%

explained_variance{iAn}.ev_bouton_pupil=ev_bouton_pupil;
explained_variance{iAn}.ev_bouton_speed=ev_bouton_speed;
explained_variance{iAn}.ev_bouton_population=ev_bouton_population;

corr_bouton_pupil=squeeze(corr_bouton_pupil(:,2,1));
corr_bouton_population=squeeze(corr_bouton_population(:,2,1));
corr_bouton_velocity=squeeze(corr_bouton_velocity(:,2,1));
corr_population_pupil=corr_population_pupil(2,1)
corr_velocity_pupil=corr_velocity_pupil(2,1);
corr_population_velocity=corr_population_velocity(2,1);
%
correlationcoef{iAn}.corr_bouton_pupil=corr_bouton_pupil;
correlationcoef{iAn}.corr_bouton_population=corr_bouton_population;
correlationcoef{iAn}.corr_bouton_velocity=corr_bouton_velocity;
correlationcoef{iAn}.corr_population_pupil=corr_population_pupil;
correlationcoef{iAn}.corr_velocity_pupil=corr_velocity_pupil;
correlationcoef{iAn}.corr_population_velocity=corr_population_velocity;
end
%% scatter plot
for iAn=1:23; 
    corr_pv_points(iAn)=correlationcoef{iAn}.corr_population_velocity;
    corr_pp_points(iAn)=correlationcoef{iAn}.corr_population_pupil;
end;
%%
figure,
subplot(221);
plot(corr_pv_points,corr_pp_points,'.','markersize',10,'color','k');
hold on
plot([-1 1],[-1 1],'--k');
axis tight
axis square
ylim([-0.2 1]);
xlim([-0.2 1]);
ylabel('Correlation Population-Velocity');
xlabel('Correlation Population-Pupil');

set(gca,'xtick',[-1:0.2:1],'ytick',[-1:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences\']; 
%print(gcf,'-dpdf',[filepathanalysis, '\fig_corrcoef_population_velocity_pupil_revision.pdf']);

%% find the same animals
%%
figure
addpath 'G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals';
clear animals_array
pattern_string='KS'
[animals_array]=find_sameanimals(newdir',pattern_string)

% [x y]=find(animals_array==iAn)
% animals_array(iAn,~isnan(animals_array(iAn,:)))

for iAn=1:23
    [x y]=find(animals_array==iAn);
    x_animal(iAn)=x
    y_session(iAn)=y
end

for n_animals=1:size(animals_array,1)
corr_pv_points_sameanimals(n_animals,:)=...
    nanmean(corr_pv_points(find(x_animal(:)==n_animals)'),2);
end

for n_animals=1:size(animals_array,1)
corr_pp_points_sameanimals(n_animals,:)=...
    nanmean(corr_pp_points(find(x_animal(:)==n_animals)'),2);
end
% plot the same animals
subplot(222);
plot(corr_pv_points_sameanimals,corr_pp_points_sameanimals,'.','markersize',10,'color','k');
hold on
plot([-1 1],[-1 1],'--k');
axis tight
axis square
ylim([-0.2 1]);
xlim([-0.2 1]);
ylabel('Correlation Population-Velocity');
xlabel('Correlation Population-Pupil');

set(gca,'xtick',[-1:0.2:1],'ytick',[-1:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences\']; 
%print(gcf,'-dpdf',[filepathanalysis, '\fig_corrcoef_population_velocity_pupil_revision_sameanimals.pdf']);
%%
for n_animals=1:size(animals_array,1)
nsession{n_animals}=find(x_animal(:)==n_animals)';
end
%% scatter plot
clear dir_pref_example
clear ori_pref_example
clear ori_index_example
clear dir_index_example
clear dir_pref_example_norm
clear ori_pref_example_norm
clear ori_index_example_norm
clear dir_index_example_norm

dir_pref_example=[stats_boutons{iAn}.dir_vector_sum_angle_degree];
ori_pref_example=[stats_boutons{iAn}.ori_vector_sum_angle_degree];
ori_index_example=[stats_boutons{iAn}.ori_vector_sum_tune];
dir_index_example=[stats_boutons{iAn}.dir_vector_sum_tune];

dir_pref_example_norm=[stats_boutons_normalized{iAn}.dir_vector_sum_angle_degree];
ori_pref_example_norm=[stats_boutons_normalized{iAn}.ori_vector_sum_angle_degree];
ori_index_example_norm=[stats_boutons_normalized{iAn}.ori_vector_sum_tune]
dir_index_example_norm=[stats_boutons_normalized{iAn}.dir_vector_sum_tune]

%% ALL EXPERIMENTS scatter plot

clear corr_bouton_pupil_total
clear corr_bouton_velocity_total
clear corr_population_pupil_total
clear corr_velocity_pupil_total
clear corr_population_velocity_total
clear corr_bouton_population_total

for iAn=1:23
corr_bouton_pupil_total{iAn}=[correlationcoef{iAn}.corr_bouton_pupil];
corr_bouton_velocity_total{iAn}=[correlationcoef{iAn}.corr_bouton_velocity];
corr_population_pupil_total{iAn}=[correlationcoef{iAn}.corr_population_pupil];
corr_velocity_pupil_total{iAn}=[correlationcoef{iAn}.corr_velocity_pupil];
corr_population_velocity_total{iAn}=[correlationcoef{iAn}.corr_population_velocity];
corr_bouton_population_total{iAn}=[correlationcoef{iAn}.corr_bouton_population];

direction_preference{iAn}=[stats_boutons{iAn}.dir_vector_sum_angle_degree];
orientation_preference{iAn}=[stats_boutons{iAn}.ori_vector_sum_angle_degree];
orientation_index{iAn}=[stats_boutons{iAn}.ori_vector_sum_tune];
direction_index{iAn}=[stats_boutons{iAn}.dir_vector_sum_tune];

end

%% TRANSFORM TO ONE ANIMAL
addpath 'G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals';
clear animals_array
pattern_string='KS'
[animals_array]=find_sameanimals(newdir',pattern_string)

for iAn=1:23
    [x y]=find(animals_array==iAn);
    x_animal(iAn)=x
    y_session(iAn)=y
end
%
for n_animals=1:size(animals_array,1)
    same_animal_sessions=find(x_animal(:)==n_animals)';
    clear corr_bouton_velo_temporary
    clear corr_bouton_pupil_temporary
    clear corr_bouton_population_temporary
    clear direction_index_temporary
    clear direction_preference_temporary
    j=1;
    
    while j<=length(same_animal_sessions)        
    corr_bouton_velo_temporary{j}=corr_bouton_velocity_total{same_animal_sessions(j)};
    corr_bouton_pupil_temporary{j}=corr_bouton_pupil_total{same_animal_sessions(j)};
    corr_bouton_population_temporary{j}=corr_bouton_population_total{same_animal_sessions(j)};
    direction_preference_temporary{j}=direction_preference{same_animal_sessions(j)}';
    direction_index_temporary{j}=direction_index{same_animal_sessions(j)}';
    j=j+1
    end
    
    corr_bouton_velocity_sameanimals{n_animals}=cell2mat(corr_bouton_velo_temporary');
    corr_bouton_pupil_sameanimals{n_animals}=cell2mat(corr_bouton_pupil_temporary');
    corr_bouton_population_sameanimals{n_animals}=cell2mat(corr_bouton_population_temporary');
    direction_preference_sameanimals{n_animals}=cell2mat(direction_preference_temporary');
    direction_index_sameanimals{n_animals}=cell2mat(direction_index_temporary');

end
%% PLOT INDIVIDUAL ANIMALS

figure
clf
subplot(224)
ed_pref = 0:22.5:360;
clear m_pref_dir
for n_animals=1:size(animals_array,1)
    clear dir_pref_example2
    clear corr_bouton_pupil2
    corr_bouton_pupil2=corr_bouton_pupil_sameanimals{n_animals}
    dir_pref_example2=direction_preference_sameanimals{n_animals}
    
list=find(dir_pref_example2>0)
clear ed
clear y_mu
[ed, y_mu, y_s] = binSamples(dir_pref_example2(list), corr_bouton_pupil2(list), ed_pref);
ed_dir=ed;
%plot(dir_pref_example2(:), corr_bouton_pupil2(:),'.','color',[0.5 0.5 0.5]);
%hold on
plot(ed_dir,y_mu,'color','k','linewidth',1);
hold on

m_pref_dir(n_animals,:)=y_mu;
ed_dir=ed;
end
subplot(224);
hold on
plot(ed_dir,nanmean(m_pref_dir,1),'color','r','linewidth',2);
set(gca,'xtick',ed_dir(1:4:end),'xticklabel',ed_pref(1:4:end),'ytick',[-0.6:0.2:1],'xticklabel',ed_pref(1:4:end),'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4);      
axis tight;
xlabel('DSI Pref Awake(n=10)');
ylabel('Correlation Bouton-Pupil');
ylim([-0.6 1]);
axis square

subplot(223)
ed_pref = 0:22.5:360;
clear m_pref_dir
for n_animals=1:size(animals_array,1)
    clear dir_pref_example2
    clear corr_bouton_pupil2
    corr_bouton_population2=corr_bouton_population_sameanimals{n_animals}
    dir_pref_example2=direction_preference_sameanimals{n_animals}
    
list=find(dir_pref_example2>0)
clear ed
clear y_mu
[ed, y_mu, y_s] = binSamples(dir_pref_example2(list), corr_bouton_population2(list), ed_pref);
ed_dir=ed;
%plot(dir_pref_example2(:), corr_bouton_pupil2(:),'.','color',[0.5 0.5 0.5]);
%hold on
plot(ed_dir,y_mu,'color','k','linewidth',1);
hold on

m_pref_dir(n_animals,:)=y_mu;
ed_dir=ed;
end
subplot(223);
hold on
plot(ed_dir,nanmean(m_pref_dir,1),'color','r','linewidth',2);
set(gca,'xtick',ed_dir(1:4:end),'xticklabel',ed_pref(1:4:end),'ytick',[-0.6:0.2:1],'xticklabel',ed_pref(1:4:end),'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4);      
axis tight;
xlabel('DSI Pref Awake(n=10)');
ylabel('Correlation Bouton-Population');
ylim([-0.6 1]);
axis square

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences\']; 
%print(gcf,'-dpdf',[filepathanalysis, '\fig_corrcoef_pref_dir_corr_sameanimals.pdf']);
