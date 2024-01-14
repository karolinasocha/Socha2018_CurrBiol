%
clear all
newdirs_anesthesia={'160707_KS164_2P_KS\run03_ori12_V1_anesthesia',...
    '160712_KS167_2P_KS\run03_ori12_V1_anesthesia',...
    '160621_KS166_2P_KS\run03_ori12_V1_anesthetized2',...
    '170110_KS173_2P_KS\run03_ori12_V1_anesthesia',...
    '170106_KS174_2P_KS\run03_ori12_V1_anesthesia',...
    '170110_KS174_2P_KS\run03_ori12_V1_anesthesia'};
%%

analys_path='G:\mousebox\analysis\';

for iAn=1:size(newdirs_anesthesia,2)
    expt_anesthesia{iAn} = frGetExpt(newdirs_anesthesia{iAn});
    expt2_anesthesia{iAn} =doLoadStimLogs3(expt_anesthesia{iAn});
    epochs{iAn} = expt2_anesthesia{iAn}.frames.epochs;
    stims{iAn} = expt2_anesthesia{iAn}.frames.stims;
    tcss_anesthesia{iAn}=load([analys_path,newdirs_anesthesia{iAn},'\tcs_handseg.mat']);
    tcss_anesthesia{iAn}=tcss_anesthesia{iAn}.tcs_handseg;
end

%%
clear max_resp

for iAn=1:size(newdirs_anesthesia,2)
    clear trials
    clear tc
tc=zscore(tcss_anesthesia{iAn}.ratio_vis);
stims2=expt2_anesthesia{iAn}.frames.stims;
%number_all_boutons(iAn)=size(tcss_anesthesia{iAn}.ratio,2);
%number_responsive_boutons(iAn)=size(tcss_anesthesia{iAn}.ratio_vis,2);
[av trials]=tcEpochAverage2(tc,stims2);

for iB=1:size(trials,2)
for itrials=1:size(trials,3)
max_resp_anesthesia{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,1:end-1))));
max_blank_anesthesia{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,end))));
end
end
%trials_resp_anesthesia{iAn}(:,:,:)=squeeze(nanmean(trials,1));
end
%%
for iAn=1:6
av_boutons(iAn)=nanmean(nanmean(max_blank_anesthesia{iAn},2));
end
mean_anesthesia_peak_blank=nanmean(av_boutons);
sem_anesthesia_peak_blank=nanstd(av_boutons)/sqrt(5);

for iAn=1:6
av_stim_boutons(iAn)=nanmean(nanmean(max_resp_anesthesia{iAn},2));
end
mean_anesthesia_peak_stim=nanmean(av_stim_boutons);
sem_anesthesia_peak_stim=nanstd(av_stim_boutons)/sqrt(5);
%%
clear max_resp2_anesthesia
clear max_blank2_anesthesia
for iAn=1:size(newdirs_anesthesia,2)  
max_resp2_anesthesia{iAn}=max_resp_anesthesia{iAn}(:);
max_blank2_anesthesia{iAn}=max_blank_anesthesia{iAn}(:);
end
tt1=cell2mat(max_resp2_anesthesia');
tt1_blank=cell2mat(max_blank2_anesthesia');

figure
subplot(221);
cd2=cdfplot(tt1);
set(cd2,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_anesthesia,2)
    cd1=cdfplot(max_resp2_anesthesia{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 10])
ylim([0 1])
set(gca,'xtick',[0:2:10],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Anesthesia')


subplot(223);
cd3=cdfplot(tt1_blank);
set(cd3,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_anesthesia,2)
    cd1=cdfplot(max_blank2_anesthesia{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 10])
ylim([0 1])
set(gca,'xtick',[0:2:10],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Anesthesia Baseline')
%%

newdirs_awake={'160707_KS164_2P_KS\run03_ori12_V1_awake',...
    '160712_KS167_2P_KS\run03_ori12_V1_awake',...
    '160621_KS166_2P_KS\run03_ori12_V1_awake',...
    '170110_KS173_2P_KS\run03_ori12_V1_awake',...
    '170106_KS174_2P_KS\run03_ori12_V1',...
    '170110_KS174_2P_KS\run03_ori12_V1_awake'};

analys_path='G:\mousebox\analysis\';
%%
for iAn=1:size(newdirs_awake,2)
    expt_awake{iAn} = frGetExpt(newdirs_awake{iAn});
    expt2_awake{iAn} =doLoadStimLogs3(expt_awake{iAn});
    epochs{iAn} = expt2_awake{iAn}.frames.epochs;
    stims{iAn} = expt2_awake{iAn}.frames.stims;
    tcss_awake{iAn}=load([analys_path,newdirs_awake{iAn},'\tcs_handseg.mat']);
    tcss_awake{iAn}=tcss_awake{iAn}.tcs_handseg;
end

%%
clear max_resp_awake

for iAn=1:size(newdirs_awake,2)
    clear trials
    clear tc
    clear av
tc=zscore(tcss_awake{iAn}.ratio_vis);
stims2=expt2_awake{iAn}.frames.stims;
%number_all_boutons(iAn)=size(tcss_anesthesia{iAn}.ratio,2);
%number_responsive_boutons(iAn)=size(tcss_anesthesia{iAn}.ratio_vis,2);
[av trials]=tcEpochAverage2(tc,stims2);

for iB=1:size(trials,2)
for itrials=1:size(trials,3)
max_resp_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,1:end-1))));
max_blank_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,end))));

end
end
%trials_resp_anesthesia{iAn}(:,:,:)=squeeze(nanmean(trials,1));
end
%%
clear max_resp2_awake
clear tt1

for iAn=1:size(newdirs_awake,2)
max_resp2_awake{iAn}=max_resp_awake{iAn}(:);
max_blank2_awake{iAn}=max_blank_awake{iAn}(:);
end
tt2=cell2mat(max_resp2_awake');
tt2_blank=cell2mat(max_blank2_awake');

subplot(222);
cd2=cdfplot(tt2);
set(cd2,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd1=cdfplot(max_resp2_awake{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 10])
ylim([0 1])
set(gca,'xtick',[0:2:10],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Awake')

subplot(224);
cd3=cdfplot(tt2_blank);
set(cd3,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd1=cdfplot(max_blank2_awake{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 10])
ylim([0 1])
set(gca,'xtick',[0:2:10],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Awake Baseline')

%%
set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Socha2018_revision\anesthesia_awake_dF_comparison\']; 
%print(gcf,'-dpdf',[filepathanalysis, '\peak_response_zscore_trials_anesthesia_6expt.pdf']);
%%

clear mean_resp_awake
clear mean_blank_awake

for iAn=1:size(newdirs_awake,2)
    clear trials
    clear tc
    clear av
tc=tcss_awake{iAn}.ratio_vis;
stims2=expt2_awake{iAn}.frames.stims;
[av trials]=tcEpochAverage2(tc,stims2);

for iB=1:size(trials,2)
for itrials=1:size(trials,3)
mean_resp_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmean(trials(:,iB,itrials,1:end-1))));
mean_blank_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmean(trials(:,iB,itrials,end))));

mean_resp_awake2{iAn}(iB,itrials)=nanmean(squeeze(nanmean(trials(:,iB,itrials,1:end-1))));
mean_blank_awake2{iAn}(iB,itrials)=nanmean(squeeze(nanmean(trials(:,iB,itrials,end))));

end
end
end
%%
for iAn=1:6
    MEAN_RESP_AWAKE{iAn}=nanmean(mean_resp_awake2{iAn},2);
    MEAN_BLANK_AWAKE{iAn}=nanmean(mean_blank_awake2{iAn},2);
end
cell_mean_resp_awake=cell2mat(MEAN_RESP_AWAKE');
cell_mean_blank_awake=cell2mat(MEAN_BLANK_AWAKE');
%%
%% 
clear mean_resp_anesthesia
clear mean_blank_anesthesia

for iAn=1:size(newdirs_anesthesia,2)
    clear trials
    clear tc
    clear av
tc=tcss_anesthesia{iAn}.ratio_vis;
stims2=expt2_anesthesia{iAn}.frames.stims;
[av trials]=tcEpochAverage2(tc,stims2);

for iB=1:size(trials,2)
for itrials=1:size(trials,3)
mean_resp_anesthesia{iAn}(iB,itrials)=nanmean(squeeze(nanmean(trials(:,iB,itrials,1:end-1))));
mean_blank_anesthesia{iAn}(iB,itrials)=nanmean(squeeze(nanmean(trials(:,iB,itrials,end))));

max_resp_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,1:end-1))));
max_blank_awake{iAn}(iB,itrials)=nanmax(squeeze(nanmax(trials(:,iB,itrials,end))));

end
end
end
%%
for iAn=1:6
    MEAN_RESP_ANSTH{iAn}=nanmean(mean_resp_anesthesia{iAn},2);
    MEAN_BLANK_ANSTH{iAn}=nanmean(mean_blank_anesthesia{iAn},2);
end
cell_mean_resp_ansth=cell2mat(MEAN_RESP_ANSTH');
cell_mean_blank_ansth=cell2mat(MEAN_BLANK_ANSTH');
%%

%% PLOT DF/F
figure
clear mean_resp2_awake
clear tt1
% AWAKE
for iAn=1:size(newdirs_awake,2)
mean_resp2_awake{iAn}=mean_resp_awake{iAn}(:);
mean_blank2_awake{iAn}=mean_blank_awake{iAn}(:);

mean_mean_resp_awake{iAn}=nanmean(mean_resp_awake{iAn},2);
mean_mean_blank_awake{iAn}=nanmean(mean_blank_awake{iAn},2);

end
%%
M_resp_awake=cell2mat(mean_mean_resp_awake');
nanstd(M_resp_awake)/sqrt(5)

M_blank_awake=cell2mat(mean_mean_blank_awake');

nanstd(M_blank_awake)/sqrt(5)
%%
clear tt2
clear tt2_blank
tt2_awake=cell2mat(mean_resp2_awake');
tt2_blank_awake=cell2mat(mean_blank2_awake');

subplot(221);
cd2=cdfplot(tt2_awake);
set(cd2,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd1=cdfplot(mean_resp2_awake{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 400])
ylim([0 1])
set(gca,'xtick',[0:50:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Awake');

subplot(223);
cd4=cdfplot(tt2_blank_awake);
set(cd4,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd5=cdfplot(mean_blank2_awake{iAn});
    set(cd5,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 400])
ylim([0 1])
set(gca,'xtick',[0:50:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Awake Baseline');


% ANESTHESIA
clear tt2
clear tt2_blank
for iAn=1:size(newdirs_awake,2)
mean_resp2_anesthesia{iAn}=mean_resp_anesthesia{iAn}(:);
mean_blank2_anesthesia{iAn}=mean_blank_anesthesia{iAn}(:);

mean_mean_resp_ansth{iAn}=nanmean(mean_resp_anesthesia{iAn},2);
mean_mean_blank_ansth{iAn}=nanmean(mean_blank_anesthesia{iAn},2);

end

M_ansth=cell2mat(mean_mean_resp_ansth');
nanstd(M_ansth)/sqrt(5)
M_blank_ansth=cell2mat(mean_mean_blank_ansth');
nanstd(M_blank_ansth)/sqrt(5)
%%
[ranksumpval,~,ranksumstat_tmp]=ranksum(rmmissing(M_blank_ansth(:)),rmmissing(M_blank_awake(:)));
ranksum(rmmissing(M_blank_ansth(:)),rmmissing(M_blank_awake(:)))

%%

[ranksumpval,~,ranksumstat_tmp]=ranksum(rmmissing(M_ansth(:)),rmmissing(M_awake(:)));

[ranksumpval,~,ranksumstat_tmp]=ranksum(rmmissing(M_blank_ansth(:)),rmmissing(M_blank(:)))    
ranksum(M_blank(:),M_blank_ansth(:))
[ranksumpval,~,ranksumstat_tmp]=ttest2(rmmissing(M_blank_ansth(:)),rmmissing(M_blank(:)))    

%%
tt2_ansth=cell2mat(mean_resp2_anesthesia');
tt2_ansth_blank=cell2mat(mean_blank2_anesthesia');


subplot(222);
cd2=cdfplot(tt2_ansth);
set(cd2,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd1=cdfplot(mean_resp2_anesthesia{iAn});
    set(cd1,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 400])
ylim([0 1])
set(gca,'xtick',[0:50:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Anesthesia');


subplot(224);
cd4=cdfplot(tt2_ansth_blank);
set(cd4,'color','r','linewidth',2);

hold on
for iAn=1:size(newdirs_awake,2)
    cd5=cdfplot(mean_blank2_anesthesia{iAn});
    set(cd5,'color',[0.5 0.5 0.5]);
end
axis tight
axis square
xlim([0 400])
ylim([0 1])
set(gca,'xtick',[0:50:400],'ytick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);       
title('Anesthesia Baseline');

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'figure6SM_data_cdf_dFF_mean_ansth_awake.pdf']);
%%

 [ranksumpval,~,ranksumstat_tmp]=ranksum(rmmissing(tt2_ansth_blank(:)),rmmissing(tt2_blank_awake(:)));

    
    