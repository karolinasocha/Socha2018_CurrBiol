%% PLOT FIGURE 4; qq plots of fraction of experiments
newdir={'160707_KS164_2P_KS\run03_ori12_V1_awake',...
    '160712_KS167_2P_KS\run03_ori12_V1_awake',...
    '160621_KS166_2P_KS\run03_ori12_V1_awake',...
    '170110_KS173_2P_KS\run03_ori12_V1_awake',...
    '170106_KS174_2P_KS\run03_ori12_V1',...
    '170110_KS174_2P_KS\run03_ori12_V1_awake'};

clear stats_boutons_tmp
clear stats_boutons_normalized_tmp
for ii=1:size(newdir,2)
 
    clear stats_version_mean_values
    expt{ii} = frGetExpt(newdir{ii});
    expt2{ii} =doLoadStimLogs3(expt{ii});
    epochs{ii} = expt2{ii}.frames.epochs;
    stims{ii} = expt2{ii}.frames.stims;
    filepathanalysis=[expt{ii}.dirs.analysis,'\',expt{ii}.animal,'\',expt{ii}.tif_folder];
    fileNametcs=[filepathanalysis,'\stats_version_mean_values.mat'];
    load(fileNametcs)

stats_boutons_tmp{ii}=stats_version_mean_values.stats_boutons;
stats_boutons_normalized_tmp{ii}=stats_version_mean_values.stats_boutons_normalized;

end

clear stats_boutons_awake
clear stats_boutons_normalized_awake
for iAn=1:size(newdir,2)    
    stats_boutons_awake{iAn}=stats_boutons_tmp{iAn};    
    stats_boutons_normalized_awake{iAn}=stats_boutons_normalized_tmp{iAn};
end

%% transfer into arrays
stats_boutons_awake=cell2mat(stats_boutons_awake);
stats_boutons_normalized_awake=cell2mat(stats_boutons_normalized_awake);
%% orientation and direction tuning 
ori_awake=[stats_boutons_awake.ori_vector_sum_tune];
dir_awake=[stats_boutons_awake.dir_vector_sum_tune];
ori_norm_awake=[stats_boutons_normalized_awake.ori_vector_sum_tune];
dir_norm_awake=[stats_boutons_normalized_awake.dir_vector_sum_tune];

ori_pref_awake=[stats_boutons_awake.ori_vector_sum_angle_degree];
dir_pref_awake=[stats_boutons_awake.dir_vector_sum_angle_degree];

ori_pref_norm_awake=[stats_boutons_normalized_awake.ori_vector_sum_angle_degree];
dir_pref_norm_awake=[stats_boutons_normalized_awake.dir_vector_sum_angle_degree];
%%
%ori_awake=osi_awake;
ori_awake=osi_anesthesia;

%dir_awake=dsi_awake;
dir_anesthesia=dsi_anesthesia;

ori_pref_awake=pref_ori_angle_awake;
ori_pref_anesthesia=pref_ori_angle_anesthesia;

dir_pref_awake=pref_dir_angle_awake;
dir_pref_anesthesia=pref_dir_angle_anesthesia;

%%
figure

c=0;
a1 = angle(exp(j*(ori_pref_awake/180*pi+c)))/pi*180; % ANESTHESIA
a2 = angle(exp(j*(ori_pref_anesthesia/180*pi+c)))/pi*180; %AWAKE
%a1_norm = angle(exp(j*(ori_pref_norm_awake/180*pi+c)))/pi*180; % AWAKE 

v =0:10:100;
p1 = prctile(a1,v);
p2 = prctile(a2,v);
%p1_norm=prctile(a1_norm,v);

subplot(231)
plot(-200:10:200,-200:10:200,'k:');
hold on;
plot(p1,p2,'k','linewidth',2);
hold on;
%plot(p1_norm,p2,'color',[0.5 0.5 0.5],'linewidth',2);

axis square;

xlim([0 180]);
ylim([0 180]);
set(gca,'box','off','ytick',[0:45:180],'xtick',[0:45:180],'tickdir','out','fontsize',12);%,...
set(gca,'layer','top','color','none','fontsize',12,'ticklength',get(gca,'ticklength')*4);

xlabel 'Awake Raw';
ylabel 'Ansth Raw';


% ------------------------------------


c=pi/2
a1 = angle(exp(j*(dir_pref_awake/180*pi+c)))/pi*180; % AWAKE 
a2 = angle(exp(j*(dir_pref_anesthesia/180*pi+c)))/pi*180; % ANESTHESIA
%a1_norm = angle(exp(j*(dir_pref_norm_awake/180*pi+c)))/pi*180; % AWAKE 


v =0:10:100;
p1 = prctile(a1,v);
p2 = prctile(a2,v);
%p1_norm=prctile(a1_norm,v);
% plot(a1,dir_pref_awake,'.');
subplot(232)
plot(-200:10:200,-200:10:200,'k:');
hold on;
plot(p1,p2,'k','linewidth',2);
hold on;
%plot(p1_norm,p2,'color',[0.5 0.5 0.5],'linewidth',2);

axis square;

xlim([-180 180]);
ylim([-180 180]);
set(gca,'box','off','ytick',[-180:45:180],'xtick',[-180:45:180],'tickdir','out','fontsize',12);%,...
set(gca,'layer','top','color','none','fontsize',12,'ticklength',get(gca,'ticklength')*4);

xlabel 'Awake Raw';
ylabel 'Ansth Raw';

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'figure4_data_qqplots_ansth_awake_norm_5expt.pdf']);

