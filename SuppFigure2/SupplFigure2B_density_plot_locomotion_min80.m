%% this is script to generate Figure3F Stationary trials comparison
% edited 02-August-2023 by KSocha
% it has weighted probability to select session

%% plot density plots
%% load data
% in relation to pupil diameter
clear all
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
%%
%%
for iAn=1:length(newdir)
    newdir_animals_id{iAn}=newdir{iAn}(8:12);
end
animal_id_list=unique(newdir_animals_id);
animal_id=nan(length(newdir),1);
sessions_id=nan(length(newdir),1);

for iAn=1:length(animal_id_list)
    indexes=find(strcmp(newdir_animals_id, animal_id_list{iAn}));
    animal_id(indexes)=iAn;
    % session
    sessions_id(indexes)=1:length(indexes)
end

% session

%%
velocity_data=velocity;
for iAn=1:size(velocity_data,2)
        epochs{iAn}=expt2{iAn}.frames.epochs;
        stims{iAn}=expt2{iAn}.frames.stims;
        tc_boutons{iAn}=diameter_data{iAn};
        [tc_response_epochs{iAn} tc_trial_response_epochs{iAn}]=tcEpochAverage2(tc_boutons{iAn},epochs{iAn});
        [tc_response_stims{iAn} tc_trial_response_stims{iAn}]=tcEpochAverage2(tc_boutons{iAn},stims{iAn});    
        [tc_velocity_stims{iAn} tc_trial_velocity_stims{iAn}]=tcEpochAverage2(velocity{iAn},stims{iAn});
        [tc_velocity_epochs{iAn} tc_trial_velocity_epochs{iAn}]=tcEpochAverage2(velocity{iAn},epochs{iAn});
end
 
%%
clear response_norm_baseline
for iAn=1:size(tc_boutons,2)
    clear eye_tmp
    clear stims_temp
    clear epochs_temp
%eye_tmp=diam{iAn};

stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
epochs{iAn}=expt2{iAn}.frames.epochs;
stims{iAn}=expt2{iAn}.frames.stims;
epochs_temp=expt2{iAn}.frames.epochs;
%[av_eye_epochs{iAn} trials_eye_epochs{iAn}]=tcEpochAverage2(eye_tmp,epochs_temp);% time missing 

stims_temp=expt2{iAn}.frames.stims;
%[av_eye_stims{iAn} trials_eye_stims{iAn}]=tcEpochAverage2(eye_tmp,stims_temp);% time missing 

[av_velo{iAn} trials_velocity_stims{iAn}]=tcEpochAverage2(velocity{iAn},stims_temp);
clear eye_tmp
clear epochs_temp
clear stims_temp
clear stimulus_resp_tmp

end

%%
clear resp_loc_velo_trials
clear resp_still_velo_trials
clear resp_nondef_velo_trials
clear resp_still_eye_trials 
clear resp_nondef_eye_trials
clear resp_loc_eye_trials
clear loctrial
clear stilltrial

%
addpath 'G:\mousebox\code\mouselab\users\karolina\Attentional_Engagement_Locomotion'

    for iAn=1:size(velocity,2);

vtime=1:length(velocity{iAn});
vval=smoothdata(velocity{iAn},'sgolay'); % used to be smooth
clear tstart
clear tstop
% locthresh=1;
% stillthresh=0.5;
% perthresh_stationary=0.95;
% perthresh=0.05;

locthresh=1;
stillthresh=0.25;
perthresh_stationary=0.95;
perthresh=0.8;

tstart = nan(size(expt2{iAn}.frames.stims));
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

% [loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials(1:length(vval),vval,...
%     tstart,tstop,locthresh,stillthresh,perthresh);

  [loctrial{iAn},stilltrial{iAn}] = select_locomotion_trials_val_stationary(1:length(vval),vval,...
      tstart,tstop,locthresh,stillthresh,perthresh,perthresh_stationary);

[loc_trial loc_stim_tmp]=find(loctrial{iAn}(:,1:end-1)==1);
[still_trial still_stim_tmp]=find(stilltrial{iAn}(:,1:end-1)==1);
nondefined{iAn}= ~(stilltrial{iAn}(:,:)+loctrial{iAn}(:,:));
[nondef_trial nondef_stim_tmp]=find(nondefined{iAn}(:,1:end-1)==1);
% select pupil diameter eye after standarization without reordering
% stimulus

tc_in=tc_trial_response_stims{iAn};

resp_loc_resp_trials{iAn}=tc_in(:,:,loc_trial',loc_stim_tmp');
resp_still_resp_trials{iAn}=tc_in(:,:,still_trial',still_stim_tmp');
resp_nondef_resp_trials{iAn}=tc_in(:,:,nondef_trial',nondef_stim_tmp');

%
stim_deg_loc{iAn}=stimulus{iAn}(loc_stim_tmp);
stim_deg_still{iAn}=stimulus{iAn}(still_stim_tmp);
stim_deg_nondef{iAn}=stimulus{iAn}(nondef_stim_tmp);

clear vtime
clear vval
clear tstart
clear tstop
clear tc_loc_tmp
clear tc_still_tmp
clear tc_nondef_tmp
clear loc_stim_tmp
clear still_stim_tmp
clear nondef_stim_tmp
clear loc_trial
clear still_trial
clear nondef_trial
    end

    %%

clear av_nondef_eye_trials
clear av_still_eye_trials
clear av_loc_eye_trials

    for iAn = 1:size(resp_loc_resp_trials,2)
    
    av_loc_eye_trials{iAn}=squeeze(nanmean(nanmean(resp_loc_resp_trials{iAn},1),2));     
    av_still_eye_trials{iAn}=squeeze(nanmean(nanmean(resp_still_resp_trials{iAn},1),2));
    av_nondef_eye_trials{iAn}=squeeze(nanmean(nanmean(resp_nondef_resp_trials{iAn},1),2));

    end
%%

for iAn=1:length(newdir)
    newdir_animals_id{iAn}=newdir{iAn}(8:12);
end
animal_id_list=unique(newdir_animals_id);
animal_id=nan(length(newdir),1);
sessions_id=nan(length(newdir),1);

for iAn=1:length(animal_id_list)
    indexes=find(strcmp(newdir_animals_id, animal_id_list{iAn}));
    animal_id(indexes)=iAn;
    % session
    sessions_id(indexes)=1:length(indexes)
end


%% ORIGINAL DATA CALCULATING P-VALUE AND STATS ORIGINAL HERE
% this will be compared with bootrapping results in the next sections
n_pairs = 6;
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);

ddata_temporal_loc2=nan([length(newdir), n_pairs, 10]);
ddata_nasal_loc2=nan([length(newdir), n_pairs,10]);

for iAn=1:length(newdir)
       for iistim=1:n_pairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        clear nasal_trials_loc
        clear temporal_trials_loc
        
        nasal_trials_loc = find(stim_deg_loc{iAn}==nasal_istim);
        temporal_trials_loc= find(stim_deg_loc{iAn}==temporal_istim);
        
        clear tmp_temporal
        clear tmp_nasal
        
        tmp_nasal=nanmean(av_loc_eye_trials{iAn}(nasal_trials_loc,nasal_trials_loc)',1);;
        tmp_temporal=nanmean(av_loc_eye_trials{iAn}(temporal_trials_loc,temporal_trials_loc)',1);

        ddata_temporal_loc2(iAn,iistim,1:length(tmp_temporal))=tmp_temporal;
        ddata_nasal_loc2(iAn,iistim,1:length(tmp_nasal))=tmp_nasal;
        % cell2mat((stim_deg_still))'
                
       end
end

try
[p_original_signrank, ~, stats_original] = signrank(ddata_nasal_loc2(:),ddata_temporal_loc2(:),'Tail','right');
stats_original_signrank = stats_original.signedrank;
end
[p_original_ranksum, ~, stats_original_ranksum] = ranksum(rmmissing(ddata_nasal_loc2(:)),rmmissing(ddata_temporal_loc2(:)),'Tail','right');

[~, p_original_ttest2, ~] = ttest2(rmmissing(ddata_nasal_loc2(:)),rmmissing(ddata_temporal_loc2(:)),'Tail','right');

[~, p_original_kstest2, stat_original_kstest2]=kstest2(rmmissing(ddata_nasal_loc2(:)),rmmissing(ddata_temporal_loc2(:)),'Tail','larger');

fprintf('kstest2 original p-value:  %.2e\n', p_original_kstest2);
fprintf('ttest2 original p-value:  %.2e\n', p_original_ttest2);
fprintf('ranksum original p-value:  %.2e\n', p_original_ranksum);
fprintf('signrank original p-value:  %.2e\n', p_original_signrank);

% kstest2 original p-value:  9.99e-01
% ttest2 original p-value:  1.13e-14
% ranksum original p-value:  4.09e-16
% signrank original p-value:  1.82e-14
%% BOOTSTRAPPING and CALCULATING BOOTSTRAPP P-BALS and STATS

n_runs=1000
n_sampled_animals=4
n_sampled_trials=3
n_pairs = 6;
ndrawnpoints=3

stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);

x_nasal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs,ndrawnpoints]);
x_temporal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs,ndrawnpoints]);

animals_ids_unique=unique(animal_id);

tmpx_all=[];
tmpy_all=[];
    
for iruns=1:n_runs
    
    tmpx = [];
    tmpy = [];

    clear animal_subset
    clear random_indices

    x_nasal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);
    x_temporal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);

    % # chose randomly animals
    random_indices = randperm(numel(animals_ids_unique));
    animal_subset = animals_ids_unique(random_indices(1:n_sampled_animals));
    animals_rand_selections(iruns, :) = animal_subset;

       clear nasal_trials_loc
       clear temporal_trials_loc
       clear ddata_temporal_loc
       clear ddata_nasal_loc

    for iAn=1:length(animal_subset)
        
       clear chosen_session
       clear animal_session_indexes
       clear random_session
       animal_session_indexes=find(animal_id==animal_subset(iAn));
       % added weight probability
       clear weights
       clear random_session
       weights=(1/length(animal_session_indexes))*(ones(1,numel(animal_session_indexes)));

       numSamples = 1; % number of session selected
       random_session = randsample(numel(animal_session_indexes), numSamples, true, weights);

       %random_session=randperm(numel(animal_session_indexes)); % this was
       %before change 202-08-2023

       chosen_session=animal_session_indexes(random_session(1));
%        % sessions selection
%        animal_session_indexes=find(animal_id==animal_subset(iAn));
%        random_session=randperm(numel(animal_session_indexes));
%        chosen_session=animal_session_indexes(random_session(1));
       
       for iistim=1:n_pairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_loc{iAn} = find(stim_deg_loc{chosen_session}==nasal_istim);
        temporal_trials_loc{iAn}= find(stim_deg_loc{chosen_session}==temporal_istim);
        
        ddata_temporal_loc{iAn}=nanmean(av_loc_eye_trials{chosen_session}(temporal_trials_loc{iAn},temporal_trials_loc{iAn})',1);
        ddata_nasal_loc{iAn}=nanmean(av_loc_eye_trials{chosen_session}(nasal_trials_loc{iAn},nasal_trials_loc{iAn})',1);
        % cell2mat((stim_deg_still))'
        
        tt1(iAn)=size(ddata_temporal_loc{iAn},2);
        tt2(iAn)=size(ddata_nasal_loc{iAn},2);
%         disp(sprintf('There are %d - %d and %d - %d trials.',...
%             tt1(iAn),nasal_istim,tt2(iAn),temporal_istim));
        if tt2(iAn)>0 & tt1(iAn)>0
            clear min_tr
            clear min_tr_still
            min_tr=min([tt1(iAn) tt2(iAn)]);
            %min_tr_still(iAn)=min_tr;
            x =  [];
            y =  [];
            for iresample = 1:n_sampled_trials
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_loc{iAn};
                temporaldirrawdatadistribution=ddata_temporal_loc{iAn};
                %x_nasal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs]);
                %x_temporal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs]);
                tmp_x_nasal=datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                tmp_x_temporal=datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
                
                x_nasal(iruns,iAn,iresample,iistim,:) =  tmp_x_nasal;
                x_temporal(iruns,iAn,iresample,iistim,:) = tmp_x_temporal;

                x_nasal2(iAn,iresample,iistim,:)=tmp_x_nasal;
                x_temporal2(iAn,iresample,iistim,:)=tmp_x_temporal;
            end
        end
    end
    tmpx = [tmpx,x_nasal2];
    tmpy = [tmpy,x_temporal2];
    
    tmpx_all=[tmpx_all,x_nasal2];
    tmpy_all=[tmpy_all,x_temporal2];
       
    end
    %[pval_signrank(iruns), h0_signrank(iruns)]=signrank(tmpx(:),tmpy(:),'Tail','right');
    [~, pval_bootstrap_ttest2(iruns)]=ttest2(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');
    
    try
    [pval_bootstrap, ~, stats_bootstrap] = signrank(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');
    stats_bootstraps_signrank(iruns) = stats_bootstrap.signedrank;
    pval_bootstrap_signrank(iruns)=pval_bootstrap;
    end
    try
    [pval_bootstrap_ranksum(iruns), ~, stats_bootstrap_ranksum] = ranksum(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');

    [~, pval_bootstrap_kstest2(iruns), stats_bootstrap_kstest2]=kstest2(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','larger');
    end
    
end
%%
prctile_tresh=50;
fprintf('kstest2 bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_kstest2,prctile_tresh));
fprintf('ttest2 bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ttest2,prctile_tresh));
fprintf('ranksum bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ranksum,prctile_tresh));
fprintf('signrank bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_signrank,prctile_tresh));

% fprintf('kstest2 bootstrapping 97.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_kstest2,prctile_tresh));
% fprintf('ttest2 bootstrapping 97.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ttest2,prctile_tresh));
% fprintf('ranksum bootstrapping 97.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ranksum,prctile_tresh));
% fprintf('signrank bootstrapping 97.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_signrank,prctile_tresh));

% kstest2 bootstrapping 97.5 prctile p-value:  8.67e-01
% ttest2 bootstrapping 97.5 prctile p-value:  1.00e+00
% ranksum bootstrapping 97.5 prctile p-value:  1.00e+00
% signrank bootstrapping 97.5 prctile p-value:  1.00e+00

prctile_tresh=2.5;
fprintf('kstest2 bootstrapping 2.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_kstest2,prctile_tresh));
fprintf('ttest2 bootstrapping 2.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ttest2,prctile_tresh));
fprintf('ranksum bootstrapping 2.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ranksum,prctile_tresh));
fprintf('signrank bootstrapping 2.5 prctile p-value:  %.2e\n', prctile(pval_bootstrap_signrank,prctile_tresh));

% kstest2 bootstrapping 2.5 prctile p-value:  6.20e-23
% ttest2 bootstrapping 2.5 prctile p-value:  2.06e-14
% ranksum bootstrapping 2.5 prctile p-value:  1.29e-14
% signrank bootstrapping 2.5 prctile p-value:  3.85e-30

%%
figure
vXEdge = linspace(0,100,60);
vYEdge = linspace(0,100,60);

%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(100*tmpx_all(:)) & isnan(100*tmpy_all(:)))
mHist2d = hist2d([100*tmpx_all(:),100*tmpy_all(:)],vYEdge,vXEdge)/npoints;

subplot(221)
imagesc(vXEdge,vYEdge, mHist2d)
hold all
set(gca,'YDir','normal')
%p = scatter(tmpy(:),tmpx(:),10,[1,1,1],'filled')
%alpha(p,0.3)
axis square
colorbar()
colormap(gray)
plot([0 150],[0,150],'w-')
%title(sprintf('Resampled %d times, each sample %d trials',nresampling,ndrawnpoints))
ylabel('delta pupil Nasal directions (%)')
xlabel('delta pupil Temporal directions (%)')
title('Locomotion Trials')
caxis([0 0.0025])

set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
%%print(gcf,'-dpdf',[filepathanalysis, 'FigureSuppl2Bcorrect_resampling_plot_density_stationary__bootstrapped.pdf']);
















%% PLOT HEAT DENSITY MAP

figure
vXEdge = linspace(0,150,60);
vYEdge = linspace(0,150,60);

%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(tmpx_all(:)) & isnan(tmpy_all(:)))
mHist2d = hist2d([tmpx_all(:),tmpy_all(:)],vYEdge,vXEdge)/npoints;

subplot(221)
imagesc(vXEdge,vYEdge, mHist2d)
hold all
set(gca,'YDir','normal')
%p = scatter(tmpy(:),tmpx(:),10,[1,1,1],'filled')
%alpha(p,0.3)
axis square
colorbar()
colormap(gray)
plot([0 150],[0,150],'w-')
%title(sprintf('Resampled %d times, each sample %d trials',nresampling,ndrawnpoints))
ylabel('dF Nasal directions (%)')
xlabel('dF Temporal directions (%)')
title('Locomotion Trials')
caxis([0 0.0025])

set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
print(gcf,'-dpdf',[filepathanalysis, 'Figure3Fcorrect_resampling_plot_density_locomotion_dF_bootstrapped.pdf']);


%% ORIGINAL DATA CALCULATING P-VALUE AND STATS ORIGINAL HERE
% this will be compared with bootrapping results in the next sections
n_pairs = 6;

stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);

ddata_temporal_still2=nan([length(newdir), n_pairs, 10]);
ddata_nasal_still2=nan([length(newdir), n_pairs,10]);

for iAn=1:length(newdir)
       for iistim=1:n_pairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        clear nasal_trials_still
        clear temporal_trials_still
        
        nasal_trials_still = find(stim_deg_still{iAn}==nasal_istim);
        temporal_trials_still= find(stim_deg_still{iAn}==temporal_istim);
        
        clear tmp_temporal
        clear tmp_nasal
        
        tmp_nasal=nanmean(av_still_eye_trials{iAn}(nasal_trials_still,nasal_trials_still)',1);;
        tmp_temporal=nanmean(av_still_eye_trials{iAn}(temporal_trials_still,temporal_trials_still)',1);

        ddata_temporal_still2(iAn,iistim,1:length(tmp_temporal))=tmp_temporal;
        ddata_nasal_still2(iAn,iistim,1:length(tmp_nasal))=tmp_nasal;
        % cell2mat((stim_deg_still))'
                
       end
end

[p_original_signrank, ~, stats_original] = signrank(ddata_nasal_still2(:),ddata_temporal_still2(:),'Tail','right');
stats_original_signrank = stats_original.signedrank;

[p_original_ranksum, ~, stats_original_ranksum] = ranksum(rmmissing(ddata_nasal_still2(:)),rmmissing(ddata_temporal_still2(:)),'Tail','right');

[~, p_original_ttest2, ~] = ttest2(rmmissing(ddata_nasal_still2(:)),rmmissing(ddata_temporal_still2(:)),'Tail','right');


%[~, p_original_kstest2]=kstest2(ddata_nasal_still(:),ddata_temporal_still(:),'alpha',0.05,'Tail','larger');

[~, p_original_kstest2, stat_original_kstest2]=kstest2(rmmissing(ddata_nasal_still2(:)),rmmissing(ddata_temporal_still2(:)),'Tail','larger');
fprintf('kstest2 original p-value:  %.2e\n', p_original_kstest2);
fprintf('ttest2 original p-value:  %.2e\n', p_original_ttest2);
fprintf('ranksum original p-value:  %.2e\n', p_original_ranksum);
fprintf('signrank original p-value:  %.2e\n', p_original_signrank);

% kstest2 original p-value:  9.98e-01
% ttest2 original p-value:  1.75e-29
% ranksum original p-value:  5.65e-32
% signrank original p-value:  3.05e-49

%% BOOTSTRAPPING and CALCULATING BOOTSTRAPP P-BALS and STATS

n_runs=1000
n_sampled_animals=4
n_sampled_trials=6
n_pairs = 6;
ndrawnpoints=6

tmpx_all=[];
tmpy_all=[];
    
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);

x_nasal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs,ndrawnpoints]);
x_temporal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs,ndrawnpoints]);

x_nasal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);
x_temporal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);

animals_ids_unique=unique(animal_id);

for iruns=1:n_runs
    tmpx = [];
    tmpy = [];

    clear animals_rand_selections
    clear animal_subset
    clear random_indices

    x_nasal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);
    x_temporal2=nan([n_sampled_animals,n_sampled_trials,n_pairs,ndrawnpoints]);

    % # chose randomly animals
    random_indices = randperm(numel(animals_ids_unique));
    animal_subset = animals_ids_unique(random_indices(1:n_sampled_animals));
    animals_rand_selections(iruns, :) = animal_subset;

       clear nasal_trials_still
       clear temporal_trials_still
       clear ddata_temporal_still
       clear ddata_nasal_still
       
       clear nasal_trials_still
       clear temporal_trials_still
       clear ddata_temporal_still
       clear ddata_nasal_still
       
    for iAn=1:length(animal_subset)
        
       clear chosen_session
       clear animal_session_indexes
       clear random_session
       animal_session_indexes=find(animal_id==animal_subset(iAn));
       % added weight probability
       clear weights
       clear random_session
       weights=(1/length(animal_session_indexes))*(ones(1,numel(animal_session_indexes)));

       numSamples = 1; % number of session selected
       random_session = randsample(numel(animal_session_indexes), numSamples, true, weights);

       %random_session=randperm(numel(animal_session_indexes)); % this was
       %before change 202-08-2023

       chosen_session=animal_session_indexes(random_session(1));
       % below section was before change 02-08-2023
%        animal_session_indexes=find(animal_id==animal_subset(iAn));
%        random_session=randperm(numel(animal_session_indexes));
%        chosen_session=animal_session_indexes(random_session(1));
       
       for iistim=1:n_pairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_still{iAn} = find(stim_deg_still{chosen_session}==nasal_istim);
        temporal_trials_still{iAn}= find(stim_deg_still{chosen_session}==temporal_istim);
        
        ddata_temporal_still{iAn}=nanmean(av_still_eye_trials{chosen_session}(temporal_trials_still{iAn},temporal_trials_still{iAn})',1);
        ddata_nasal_still{iAn}=nanmean(av_still_eye_trials{chosen_session}(nasal_trials_still{iAn},nasal_trials_still{iAn})',1);
        % cell2mat((stim_deg_still))'
        
        tt1(iAn)=size(ddata_temporal_still{iAn},2);
        tt2(iAn)=size(ddata_nasal_still{iAn},2);
%         disp(sprintf('There are %d - %d and %d - %d trials.',...
%             tt1(iAn),nasal_istim,tt2(iAn),temporal_istim));
        if tt2(iAn)>0 & tt1(iAn)>0
            clear min_tr
            clear min_tr_still
            min_tr=min([tt1(iAn) tt2(iAn)]);
            %min_tr_still(iAn)=min_tr;
            x =  [];
            y =  [];
            for iresample = 1:n_sampled_trials
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_still{iAn};
                temporaldirrawdatadistribution=ddata_temporal_still{iAn};
                %x_nasal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs]);
                %x_temporal = nan([n_runs, n_sampled_animals, n_sampled_trials,n_pairs]);
                tmp_x_nasal=datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                tmp_x_temporal=datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
                
                x_nasal(iruns,iAn,iresample,iistim,:) =  tmp_x_nasal;
                x_temporal(iruns,iAn,iresample,iistim,:) = tmp_x_temporal;

                x_nasal2(iAn,iresample,iistim,:)=tmp_x_nasal;
                x_temporal2(iAn,iresample,iistim,:)=tmp_x_temporal;
            end
        end
    end
    tmpx = [tmpx,x_nasal2];
    tmpy = [tmpy,x_temporal2];
    
    tmpx_all=[tmpx_all,x_nasal2];
    tmpy_all=[tmpy_all,x_temporal2];
    end
    %[pval_signrank(iruns), h0_signrank(iruns)]=signrank(tmpx(:),tmpy(:),'Tail','right');
    [~, pval_bootstrap_ttest2(iruns)]=ttest2(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');
    
    [pval_bootstrap, ~, stats_bootstrap] = signrank(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');
    stats_bootstraps_signrank(iruns) = stats_bootstrap.signedrank;
    pval_bootstrap_signrank(iruns)=pval_bootstrap;
    
    [pval_bootstrap_ranksum(iruns), ~, stats_bootstrap_ranksum] = ranksum(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','right');

    [~, pval_bootstrap_kstest2(iruns), stats_bootstrap_kstest2]=kstest2(rmmissing(tmpx(:)),rmmissing(tmpy(:)),'Tail','larger');
    
end

prctile_tresh=50;
fprintf('kstest2 bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_kstest2,prctile_tresh));
fprintf('ttest2 bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ttest2,prctile_tresh));
fprintf('ranksum bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_ranksum,prctile_tresh));
fprintf('signrank bootstrapping 50 prctile p-value:  %.2e\n', prctile(pval_bootstrap_signrank,prctile_tresh));


%%


%% PLOT HEAT DENSITY MAP

figure
vXEdge = linspace(0,100,60);
vYEdge = linspace(0,100,60);

%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(100*tmpx_all(:)) & isnan(100*tmpy_all(:)))
mHist2d = hist2d([100*tmpx_all(:),100*tmpy_all(:)],vYEdge,vXEdge)/npoints;

subplot(221)
imagesc(vXEdge,vYEdge, mHist2d)
hold all
set(gca,'YDir','normal')
%p = scatter(tmpy(:),tmpx(:),10,[1,1,1],'filled')
%alpha(p,0.3)
axis square
colorbar()
colormap(gray)
plot([0 150],[0,150],'w-')
%title(sprintf('Resampled %d times, each sample %d trials',nresampling,ndrawnpoints))
ylabel('delta pupil Nasal directions (%)')
xlabel('delta pupil Temporal directions (%)')
title('Stationary Trials')
caxis([0 0.0025])

set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\Figure3\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'FigureSuppl2Bcorrect_resampling_plot_density_stationary__bootstrapped.pdf']);

