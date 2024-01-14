% plot cdf all experiments
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
    epochs{iAn} = expt2{iAn}.frames.epochs;
    stims{iAn} = expt2{iAn}.frames.stims;

    velocity_data{iAn}=velocity_27expt_v2{iAn};
    newdir{iAn}=newdir_27expt_v2{iAn};
     if exist([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat'])
    tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_handseg.mat']);
    tcs{iAn}=tcs{iAn}.tcs_handseg;
    else
         tcs{iAn}=load([expt{iAn}.dirs.analrootpn,'\tcs_colormap.mat']);
         tcs{iAn}=tcs{iAn}.tcs_colormap;
     end
      
     
end

%%
for iAn=1:size(newdir,2)
    clear trials
    clear stims2
    clear tc
    clear av
tc=tcs{iAn}.ratio_vis;
stims2=expt2{iAn}.frames.stims;
number_all_boutons(iAn)=size(tcs{iAn}.ratio,2);
number_responsive_boutons(iAn)=size(tcs{iAn}.ratio_vis,2);
[av trials]=tcEpochAverage2(tc,stims2);
trials_resp{iAn}=squeeze(nanmean(trials,1));
end

clear peak_response_alldatapoint
clear peak_response_blank_anesthesia
clear peak_response
clear peak_response_blank

for iAn=1:size(newdir,2)  
peak_response{iAn} = squeeze(nanmean(trials_resp{iAn}(:,:,1:end-1),3));
peak_response_blank{iAn} = squeeze(nanmean(trials_resp{iAn}(:,:,end),3));
peak_response_alldatapoint{iAn} = trials_resp{iAn}(:,:,1:end-1);
peak_response_blank_alldatapoint{iAn} =trials_resp{iAn}(:,:,end);

end

clear peak_response2_alldatapoint
clear peak_response2_blank_alldatapoint
clear peak_response2
clear peak_response2_blank
for iAn=1:size(newdir,2);
peak_response2{iAn}=peak_response{iAn}(:);
peak_response2_blank{iAn}=peak_response_blank{iAn}(:);
peak_response2_alldatapoint{iAn}=peak_response_alldatapoint{iAn}(:);
peak_response2_blank_alldatapoint{iAn}=peak_response_blank_alldatapoint{iAn}(:);

end
%% PLOT CDF AMPLITUDE RESPONSE

figure,
subplot(221)
cd1=cdfplot(cell2mat(peak_response2'));
set(cd1,'color','b')
hold on
cd2=cdfplot(cell2mat(peak_response2_blank_alldatapoint'));
set(cd2,'color','k')

axis tight
axis square
xlim([0 400]);
ylim([0 1]);
xlabel('dF/F')
ylabel('Fraction')
axis square
set(gca,'xtick',[0:100:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       


set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\anesthesia_awake_dF_comparison\']; 
%print(gcf,'-dpdf',[filepathanalysis, '23expt_cdf_dFF_peak_awake.pdf']);

%% 
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
             '170107_KS174_2P_KS\run03_ori12_V1',...% dodane
             '170110_KS174_2P_KS\run03_ori12_V1_awake',...% dodane
             '170108_KS174_2P_KS\run03_ori12_V1_awake',... % dodane
             '170104_KS173_2P_KS\run03_ori12_V1',... % dodane  remove for heatmaps 100 frames
             '170110_KS173_2P_KS\run03_ori12_V1_awake'};% dodane
         
for iAn=1:size(newdir,2)
expt{iAn} = frGetExpt(newdir{iAn});
expt2{iAn} =doLoadStimLogs3(expt{iAn});
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
for iAn=1:size(newdir,2)
    clear trials
    clear stims2
    clear tc
    clear av
tc=tcs{iAn}.ratio_vis;
stims2=expt2{iAn}.frames.stims;
[av trials]=tcEpochAverage2(tc,stims2);
trials_resp{iAn}=squeeze(nanmean(trials,1));
end

clear peak_response_alldatapoint
clear peak_response_blank_anesthesia
clear peak_response
clear peak_response_blank

for iAn=1:size(newdir,2)  
peak_response{iAn} = squeeze(nanmean(trials_resp{iAn}(:,:,1:end-1),3));
peak_response_blank{iAn} = squeeze(nanmean(trials_resp{iAn}(:,:,end),3));
peak_response_alldatapoint{iAn} = trials_resp{iAn}(:,:,1:end-1);
peak_response_blank_alldatapoint{iAn} =trials_resp{iAn}(:,:,end);

end

clear peak_response2_alldatapoint
clear peak_response2_blank_alldatapoint
clear peak_response2
clear peak_response2_blank
for iAn=1:size(newdir,2);
peak_response2{iAn}=peak_response{iAn}(:);
peak_response2_blank{iAn}=peak_response_blank{iAn}(:);
peak_response2_alldatapoint{iAn}=peak_response_alldatapoint{iAn}(:);
peak_response2_blank_alldatapoint{iAn}=peak_response_blank_alldatapoint{iAn}(:);
end
%% PLOT CDF AMPLITUDE RESPONSE

figure,
subplot(221)
cd1=cdfplot(cell2mat(peak_response2'));
set(cd1,'color','b')
hold on
cd2=cdfplot(cell2mat(peak_response2_blank_alldatapoint'));
set(cd2,'color','k')

axis tight
axis square
xlim([0 400]);
ylim([0 1]);
xlabel('dF/F')
ylabel('Fraction')
axis square
set(gca,'xtick',[0:100:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(gca,'ticklength')*4)       


set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\anesthesia_awake_dF_comparison\']; 
%print(gcf,'-dpdf',[filepathanalysis, '37expt_cdf_dFF_peak_awake.pdf']);
%%

pathname='G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\anesthesia_awake_dF_comparison\';
%%
clear sampleA_ansth
clear sampleB_awake

sampleA_ansth=cell2mat(peak_response2_anesthesia');
sampleB_awake=cell2mat(peak_response2');
[p,h]=ranksum(sampleA_ansth,sampleB_awake)

clear sampleA_ansth
clear sampleB_awake

sampleA_ansth=cell2mat(peak_response2_awake');
sampleB_awake=cell2mat(peak_response2');

[p,h]=ranksum(sampleA_ansth,sampleB_awake)

clear sampleA_ansth
clear sampleB_awake

sampleA_ansth=cell2mat(peak_response2_anesthesia');
sampleB_awake=cell2mat(peak_response2_awake');
[p,h]=ranksum(sampleA_ansth,sampleB_awake)
