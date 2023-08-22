%% plot density plots
%% load data
% in relation to pupil diameter
clear all
filepathanalysis='G:\mousebox\code\mouselab\users\karolina\figure_quantif_behaviour_eye_velocity_selectedexpt';

%%
%load data 37 expt; 8240 boutons
% load data to remove gain EXPERIMENTS COMBINED
clear all
newdir={'140623_KS092_2P_KS\run01_ori_ds_V1',...
    '140623_KS093_2P_KS\run01_ori_ds_V1',...
    '150420_KS133_2P_KS\run01_ori_ds_V1',...
    '150915_KS145_2P_KS\run01_ori_ds_V1_full',...
    '140703_KS099_2P_KS\run01_ori_ds_V1',...
    '141209_KS126_2P_KS\run01_ori_ds_V1_reversed',... %remove for heatmaps  has 9 stims problematic
    '140810_KS103_2P_KS\run01_ori_ds_V1_full',...
    '160105_KS154_2P_KS\run01_orids_V1',... %remove for heatmaps 9 stims
    '151007_KS145_2P_KS\run012_ori_ds_V1_full',...dodane
  '140808_KS092_2P_KS\run02_ori_ds_V1_full',... % remove for heatmaps 9 stims
  '140808_KS093_2P_KS\run02_ori_ds_V1_full',...
  '150420_KS133_2P_KS\run02_ori_ds_V1',...
  '150917_KS145_2P_KS\run02_ori_ds_V1_full',...
  '141215_KS126_2P_KS\run02_ori_ds_V1_reversed',... % remove for heatmaps 9 stims
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
             '170107_KS174_2P_KS\run03_ori12_V1',... % remove for heatmaps zmienic potem na 170110
             '170110_KS173_2P_KS\run03_ori12_V1_awake',...
             '151007_KS145_2P_KS\run034_ori_ds_V1_full',... % dodane
             '151007_KS145_2P_KS\run03_ori_ds_V1_full',...% dodane
             ... %'170107_KS174_2P_KS\run03_ori12_V1',...% dodane
             '170110_KS174_2P_KS\run03_ori12_V1_awake',...% dodane
             '170108_KS174_2P_KS\run03_ori12_V1_awake',... % dodane
             '170104_KS173_2P_KS\run03_ori12_V1',... % dodane  remove for heatmaps 100 frames
             '170110_KS173_2P_KS\run03_ori12_V1_awake'};% dodane
         
for iAn=1:size(newdir,2)
expt{iAn} = frGetExpt(newdir{iAn});
expt2{iAn} =doLoadStimLogs3(expt{iAn});
stimulus{iAn}=[expt2{iAn}.prot.pars.ori];
vel{iAn} = load([expt{iAn}.dirs.analrootpn,'\','velo_final.mat']);
velocity{iAn} = vel{iAn}.velo_final(:,2);
epochs{iAn} = expt2{iAn}.frames.epochs;
stims{iAn} = expt2{iAn}.frames.stims;
     if exist([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat'])
    tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat']);
    tcs{iAn}=tcs{iAn}.tcs_handseg;
    else
         tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_colormap.mat']);
         tcs{iAn}=tcs{iAn}.tcs_colormap;
     end
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
for iAn=1:size(tcs,2)
        tc_boutons{iAn}=tcs{iAn}.ratio_vis;
        [tc_response_epochs{iAn} tc_trial_response_epochs{iAn}]=tcEpochAverage2(tc_boutons{iAn},epochs{iAn});
        [tc_response_stims{iAn} tc_trial_response_stims{iAn}]=tcEpochAverage2(tc_boutons{iAn},stims{iAn});    
        [tc_velocity_stims{iAn} tc_trial_velocity_stims{iAn}]=tcEpochAverage2(velocity{iAn},stims{iAn});
        [tc_velocity_epochs{iAn} tc_trial_velocity_epochs{iAn}]=tcEpochAverage2(velocity{iAn},epochs{iAn});
end
 
%%
clear response_norm_baseline
for iAn=1:size(tcs,2)
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
% stimulus_resp_tmp=trials_eye_stims{iAn};
% for iTrial=1:size(stimulus_resp_tmp,3);
%     for iStim=1:size(stimulus_resp_tmp,4);
%         response_norm_baseline{iAn}(:,:,iTrial,iStim)=(stimulus_resp_tmp(:,:,iTrial,iStim)-nanmedian(diameter{iAn}))./nanmedian(diameter{iAn});
%     end
% end

end
%%
%trials_eye_stims=response_norm_baseline

%%
% clear animals_array
%  addpath 'G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals'
%  pattern_string='KS'
% [animals_array]=find_sameanimals(newdir,pattern_string)
% 
% [x y]=find(animals_array==iAn)
% %animals_array(iAn,~isnan(animals_array(iAn,:)))
% 
% for iAn=1:size(newdir,2)
%     [x y]=find(animals_array==iAn);
%     x_animal(iAn)=x;
%     y_session(iAn)=y;
% end
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
locthresh=1;
stillthresh=1;
perthresh=0.95;

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
locthresh=1;
stillthresh=1;
perthresh_stationary=0.95;
perthresh=0.05;
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

%% prepare array with values depending on behaviour condition

% nasal_trials_mixed{iAn} = find(stim_deg_nondef{iAn}==nasal_istim);
% temporal_trials_mixed{iAn}= find(stim_deg_nondef{iAn}==temporal_istim);
% ddata_temporal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(temporal_trials_mixed{iAn},temporal_trials_mixed{iAn})',1);
% ddata_nasal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(nasal_trials_mixed{iAn},nasal_trials_mixed{iAn})',1);

% nasal_trials_still{iAn} = find(stim_deg_still{iAn}==nasal_istim);
% temporal_trials_still{iAn}= find(stim_deg_still{iAn}==temporal_istim);
% ddata_temporal_still{iAn}=nanmean(av_still_eye_trials{iAn}(temporal_trials_still{iAn},temporal_trials_still{iAn})',1);
% ddata_nasal_still{iAn}=nanmean(av_still_eye_trials{iAn}(nasal_trials_still{iAn},nasal_trials_still{iAn})',1);

% nasal_trials_locomotion{iAn} = find(stim_deg_loc{iAn}==nasal_istim);
% temporal_trials_locomotion{iAn}= find(stim_deg_loc{iAn}==temporal_istim);
% ddata_temporal_locomotion{iAn}=nanmean(av_loc_eye_trials{iAn}(temporal_trials_locomotion{iAn},temporal_trials_locomotion{iAn})',1);
% ddata_nasal_locomotion{iAn}=nanmean(av_loc_eye_trials{iAn}(nasal_trials_locomotion{iAn},nasal_trials_locomotion{iAn})',1);
%% create additional matrix about animal_id; sessions_id; stimulus value

for ii=1:size(stim_deg_still,2)
    animal_array{ii}=size(stim_deg_still{iAn})
end
%% example one animal ALL SESSIONS STATIONARY
clear tmpx
clear tmpy
rng('default');
rng(123);
ndrawnpoints = 3;
nresampling = 1000;
npairs = 6;
tmpx = [];
tmpy = [];
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);

for iAn = 1:length(newdir)
    clear nasal_trials_still
    clear temporal_trials_still
    clear ddata_temporal_still
    clear ddata_nasal_still
    clear x_nasal
    clear x_temporal
    x_nasal = zeros([ndrawnpoints,nresampling,npairs]);
    x_temporal = zeros([ndrawnpoints,nresampling,npairs]);
    x_nasal(:) = nan;
    x_temporal(:) = nan;
    
    for iistim=1:npairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_still{iAn} = find(stim_deg_still{iAn}==nasal_istim);
        temporal_trials_still{iAn}= find(stim_deg_still{iAn}==temporal_istim);
        ddata_temporal_still{iAn}=nanmean(av_still_eye_trials{iAn}(temporal_trials_still{iAn},temporal_trials_still{iAn})',1);
        ddata_nasal_still{iAn}=nanmean(av_still_eye_trials{iAn}(nasal_trials_still{iAn},nasal_trials_still{iAn})',1);
        % cell2mat((stim_deg_still))'
        clear tt1
        clear tt2
        
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
            for iresample = 1:nresampling
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_still{iAn};
                temporaldirrawdatadistribution=ddata_temporal_still{iAn};
                
                x_nasal(:,iresample,iistim) =  datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                x_temporal(:,iresample,iistim) = datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
            end
        end
    end
    tmpx = [tmpx,x_nasal];
    tmpy = [tmpy,x_temporal];
    
end
%%
%%
figure
% vXEdge = linspace(-0.8,0.8,75);
% vYEdge = linspace(-0.8,0.8,75);
vXEdge = linspace(0,150,60);
vYEdge = linspace(0,150,60);
%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(tmpx(:)) & isnan(tmpy(:)))
mHist2d = hist2d([tmpx(:),tmpy(:)],vYEdge,vXEdge)/npoints;

subplot(223)
imagesc(vXEdge,vYEdge, mHist2d)
hold all
set(gca,'YDir','normal')
%p = scatter(tmpy(:),tmpx(:),10,[1,1,1],'filled')
%alpha(p,0.3)
axis square
colorbar()
colormap(gray)
%plot([-0.8,0.8],[-0.8,0.8],'w-')
plot([0,150],[0,150],'w-')
ylabel('Stationary')
set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       
caxis([0 0.005])
%title(sprintf('Resampled %d times, each sample %d trials',nresampling,ndrawnpoints))

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'resampling_plot_density_allconditions_gray_colormap.pdf']);
%%
help kstest2
%[h_still, p_still]=kstest2(tmpx(:), tmpy(:));
indx=find(isnan(tmpx(:))==1);
tmpx(indx)=[];

clear indx
indx=find(isnan(tmpy(:))==1);
tmpy(indx)=[];
% tmpx - nasal
% tmpy - temporal
[p_still_signrank, h_still_signrank]=signrank(tmpx(:), tmpy(:),'Tail','right');

%% LOCOMOTION ALL SESSIONS
% 
rng('default');
rng(123);
ndrawnpoints = 3;
nresampling = 1000;
npairs = 6;
tmpx = [];
tmpy = [];
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);
for iAn = 1:length(newdir)
    clear nasal_trials_still
    clear temporal_trials_still
    clear ddata_temporal_still
    clear ddata_nasal_still
    clear x_nasal
    clear x_temporal
    x_nasal = zeros([ndrawnpoints,nresampling,npairs]);
    x_temporal = zeros([ndrawnpoints,nresampling,npairs]);
    x_nasal(:) = nan;
    x_temporal(:) = nan;
    for iistim=1:npairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_locomotion{iAn} = find(stim_deg_loc{iAn}==nasal_istim);
        temporal_trials_locomotion{iAn}= find(stim_deg_loc{iAn}==temporal_istim);
        ddata_temporal_locomotion{iAn}=nanmean(av_loc_eye_trials{iAn}(temporal_trials_locomotion{iAn},temporal_trials_locomotion{iAn})',1);
        ddata_nasal_locomotion{iAn}=nanmean(av_loc_eye_trials{iAn}(nasal_trials_locomotion{iAn},nasal_trials_locomotion{iAn})',1);

        
        
        clear tt1
        clear tt2
        
        tt1(iAn)=size(ddata_temporal_locomotion{iAn},2);
        tt2(iAn)=size(ddata_nasal_locomotion{iAn},2);
%         disp(sprintf('There are %d - %d and %d - %d trials.',...
%             tt1(iAn),nasal_istim,tt2(iAn),temporal_istim));
        if tt2(iAn)>1 & tt1(iAn)>1
            clear min_tr
            clear min_tr_still
            min_tr=min([tt1(iAn) tt2(iAn)]);
            %min_tr_still(iAn)=min_tr;
            x =  [];
            y =  [];
            for iresample = 1:nresampling
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_locomotion{iAn};
                temporaldirrawdatadistribution=ddata_temporal_locomotion{iAn};
                
                x_nasal(:,iresample,iistim) =  datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                x_temporal(:,iresample,iistim) = datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
            end
        end
    end
    tmpx = [tmpx,x_nasal];
    tmpy = [tmpy,x_temporal];
    
end
%%
clear indx_nasal
clear indx_temporal
tmpx2=tmpx(:);
tmpy2=tmpy(:);
isnan(tmpx2);
indx_nasal=find(~isnan(tmpx2(:))==0); 
indx_temporal=find(~isnan(tmpy2(:))==0); 

tmpx2([indx_nasal indx_temporal])=[];
tmpy2([indx_nasal indx_temporal])=[];

[h_locomotion, p_locomotion]=kstest2(tmpx2, tmpy2);
figure
c1=cdfplot(tmpx2(:));
hold on
c2=cdfplot(tmpy2(:));
set(c2,'color','k');

%%
figure
vXEdge = linspace(0,150,60);
vYEdge = linspace(0,150,60);

%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(tmpx(:)) & isnan(tmpy(:)))
mHist2d = hist2d([tmpx(:),tmpy(:)],vYEdge,vXEdge)/npoints;

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
ylabel('Locomotion + Mixed')
caxis([0 0.005])

set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\density_plots_dF\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'resampling_plot_density_locomotion_stationary_dF_samescale.pdf']);

%%
   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'resampling_plot_density_locomotion_differentscale.pdf']);

%% MIXED ALL SESSIONS
% 
rng('default');
rng(123);
ndrawnpoints = 3;
nresampling = 1000;
npairs = 6;
tmpx = [];
tmpy = [];
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);
for iAn = 1:length(newdir)
    clear nasal_trials_still
    clear temporal_trials_still
    clear ddata_temporal_still
    clear ddata_nasal_still
    clear x_nasal
    clear x_temporal
    x_nasal = zeros([ndrawnpoints,nresampling,npairs]);
    x_temporal = zeros([ndrawnpoints,nresampling,npairs]);
    x_nasal(:) = nan;
    x_temporal(:) = nan;
    for iistim=1:npairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_mixed{iAn} = find(stim_deg_nondef{iAn}==nasal_istim);
        temporal_trials_mixed{iAn}= find(stim_deg_nondef{iAn}==temporal_istim);
        ddata_temporal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(temporal_trials_mixed{iAn},temporal_trials_mixed{iAn})',1);
        ddata_nasal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(nasal_trials_mixed{iAn},nasal_trials_mixed{iAn})',1);

        
        clear tt1
        clear tt2
        
        tt1(iAn)=size(ddata_temporal_mixed{iAn},2);
        tt2(iAn)=size(ddata_nasal_mixed{iAn},2);
%         disp(sprintf('There are %d - %d and %d - %d trials.',...
%             tt1(iAn),nasal_istim,tt2(iAn),temporal_istim));
        if tt2(iAn)>1 & tt1(iAn)>1
            clear min_tr
            clear min_tr_still
            min_tr=min([tt1(iAn) tt2(iAn)]);
            %min_tr_still(iAn)=min_tr;
            x =  [];
            y =  [];
            for iresample = 1:nresampling
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_mixed{iAn};
                temporaldirrawdatadistribution=ddata_temporal_mixed{iAn};
                
                x_nasal(:,iresample,iistim) =  datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                x_temporal(:,iresample,iistim) = datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
            end
        end
    end
    tmpx = [tmpx,x_nasal];
    tmpy = [tmpy,x_temporal];
    
end
%%
[h_mixed, p_mixed]=kstest2(tmpx(:), tmpy(:));

%%
figure
vXEdge = linspace(0,250,60);
vYEdge = linspace(0,250,60);

%mHist2d = hist2d([x_nasal(:),x_temporal(:)],vYEdge,vXEdge)/length(x_nasal(:));
npoints = sum(isnan(tmpx(:)) & isnan(tmpy(:)))
mHist2d = hist2d([tmpx(:),tmpy(:)],vYEdge,vXEdge)/npoints;

subplot(222)
imagesc(vXEdge,vYEdge, mHist2d)
hold all
set(gca,'YDir','normal')
%p = scatter(tmpy(:),tmpx(:),10,[1,1,1],'filled')
%alpha(p,0.3)

colorbar()
colormap(bluered)
plot([0,250],[0,250],'w-');
axis square
%title(sprintf('Resampled %d times, each sample %d trials',nresampling,ndrawnpoints))
ylabel('Mixed')
set(gca,'ytick',[0:50:150],'xtick',[0:50:150],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       

   set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
   filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\qunatify_same_animals\']; 
   %print(gcf,'-dpdf',[filepathanalysis, 'resampling_plot_density_allconditions_normalized_data_gray.pdf']);
%%
%% MIXED ALL SESSIONS
% 
rng('default');
rng(123);
ndrawnpoints = 3;
nresampling = 1000;
npairs = 6;
tmpx = [];
tmpy = [];
stimulus_ordered=0:30:330;
stimulus_ordered_mod=mod(stimulus_ordered-60,360);
for iAn = 1:length(newdir)
    clear nasal_trials_still
    clear temporal_trials_still
    clear ddata_temporal_still
    clear ddata_nasal_still
    clear x_nasal
    clear x_temporal
    x_nasal = zeros([ndrawnpoints,nresampling,npairs]);
    x_temporal = zeros([ndrawnpoints,nresampling,npairs]);
    x_nasal(:) = nan;
    x_temporal(:) = nan;
    for iistim=1:npairs
        
        nasal_istim=stimulus_ordered_mod(iistim);
        temporal_istim=stimulus_ordered_mod(iistim+6);
        
        iStim_nasal=iistim;
        iStim_temporal=iistim+6;
        
        nasal_trials_mixed{iAn} = find(stim_deg_nondef{iAn}==nasal_istim);
        temporal_trials_mixed{iAn}= find(stim_deg_nondef{iAn}==temporal_istim);
        ddata_temporal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(temporal_trials_mixed{iAn},temporal_trials_mixed{iAn})',1);
        ddata_nasal_mixed{iAn}=nanmean(av_nondef_eye_trials{iAn}(nasal_trials_mixed{iAn},nasal_trials_mixed{iAn})',1);

        
        clear tt1
        clear tt2
        
        tt1(iAn)=size(ddata_temporal_mixed{iAn},2);
        tt2(iAn)=size(ddata_nasal_mixed{iAn},2);
%         disp(sprintf('There are %d - %d and %d - %d trials.',...
%             tt1(iAn),nasal_istim,tt2(iAn),temporal_istim));
        if tt2(iAn)>1 & tt1(iAn)>1
            clear min_tr
            clear min_tr_still
            min_tr=min([tt1(iAn) tt2(iAn)]);
            %min_tr_still(iAn)=min_tr;
            x =  [];
            y =  [];
            for iresample = 1:nresampling
                %nasaldirrawdatadistribution = GET
                %temporaldirrawdatadistribution = GET
                clear nasaldirrawdatadistribution
                clear temporaldirrawdatadistribution
                nasaldirrawdatadistribution=ddata_nasal_mixed{iAn};
                temporaldirrawdatadistribution=ddata_temporal_mixed{iAn};
                
                x_nasal(:,iresample,iistim) =  datasample(nasaldirrawdatadistribution,ndrawnpoints,'replace',true);
                x_temporal(:,iresample,iistim) = datasample(temporaldirrawdatadistribution,ndrawnpoints,'replace',true);
            end
        end
    end
    tmpx = [tmpx,x_nasal];
    tmpy = [tmpy,x_temporal];
    
end


%% MIXED if you split Mixed and locomotion

clear response_nasal_mixed
clear response_temporal_mixed
for iistim=1:6

    nasal_istim=mod(stimulus{1}(iistim)-60,360)
    temporal_istim=mod(stimulus{1}(iistim+6)-60,360)

    iStim_nasal=find(stimulus{1}(1:12)==nasal_istim);
    iStim_temporal=find(stimulus{1}==temporal_istim);

clear tt1
clear tt2
for n_animals=1:size(arr_mixed,2)
tt1(n_animals)=length(find(~isnan(unique(arr_mixed{n_animals}(iStim_nasal,:)))>0));
tt2(n_animals)=length(find(~isnan(unique(arr_mixed{n_animals}(iStim_temporal,:)))>0));
end


for n_animals=1:size(arr_mixed,2)
response_nasal_mixed{n_animals}=nan(6,15); 
response_temporal_mixed{n_animals}=nan(6,15); 

if tt1(n_animals)>0 & tt2(n_animals)>0
    clear min_tr
    clear min_tr_still
    min_tr=min([tt1(n_animals) tt2(n_animals)]);
    min_tr_still(n_animals)=min_tr;

    response_nasal_mixed{n_animals}(iistim,1:min_tr)=arr_mixed{n_animals}(iStim_nasal,randperm(max(tt1(n_animals)),min_tr))
    response_temporal_mixed{n_animals}(iistim,1:min_tr)=arr_mixed{n_animals}(iStim_temporal,randperm(max(tt2(n_animals)),min_tr))
    animal_id_mixed{n_animals}=ones(1,min_tr)*n_animals
    
end    
end
end

%clear response_temporal_still
%clear response_nasal_still

%%