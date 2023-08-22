%% figure 4 osi/dsi
%%
%load data 37 expt; 8240 boutons
% load data 36 experiments: 7908 because duplicated one experiment from
% KS174 but removed 21 August 2023
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
             ...%'170107_KS174_2P_KS\run03_ori12_V1',...% dodane removed duplicated
             '170110_KS174_2P_KS\run03_ori12_V1_awake',...% dodane
             '170108_KS174_2P_KS\run03_ori12_V1_awake',... % dodane
             '170104_KS173_2P_KS\run03_ori12_V1',... % dodane  remove for heatmaps 100 frames
             '170110_KS173_2P_KS\run03_ori12_V1_awake'};% dodane
for iAn=1:size(newdir,2)
expt{iAn} = frGetExpt(newdir{iAn});
expt2{iAn} =doLoadStimLogs3(expt{iAn});
nStim{iAn}=expt2{iAn}.info.nStim;
clear an_stat
an_stat=load([expt{iAn}.dirs.analrootpn,'\stats_version_mean_values.mat']);
stats_selected{iAn}=[an_stat.stats_version_mean_values.stats_boutons];
stats_selected_normalized{iAn}=[an_stat.stats_version_mean_values.stats_boutons_normalized];
end  
%%
stats_boutons_all=cell2mat(stats_selected)
stats_boutons_all_norm=cell2mat(stats_selected_normalized)

%%
prefosi_norm_all=[stats_boutons_all_norm.ori_vector_sum_angle_degree];
prefosi_all=[stats_boutons_all.ori_vector_sum_angle_degree];
prefdsi_norm_all=[stats_boutons_all_norm.dir_vector_sum_angle_degree];
prefdsi_all=[stats_boutons_all.dir_vector_sum_angle_degree];

osi_norm_all=[stats_boutons_all_norm.ori_vector_sum_tune];
osi_all=[stats_boutons_all.ori_vector_sum_tune];

dsi_norm_all=[stats_boutons_all_norm.dir_vector_sum_tune];
dsi_all=[stats_boutons_all.dir_vector_sum_tune];%% DSI scatter plot
%%
dsi_awake= dsi_all;
osi_awake=osi_all;
pref_ori_angle_awake=prefosi_all(find(osi_awake>0.2));
pref_dir_angle_awake=prefdsi_all(find(dsi_awake>0.2));
pref_dir_angle_awake_norm=prefdsi_norm_all(find(dsi_norm_all>0.2));
choose_color=[0 0 0];
%%
figure('name','awake')
clf
%awake_color=[110/255 150/255 200/255];
%anesthesia_color=[185/255 190/255 190/255];
awake_color=[0 0 0];
ax=[];
ax(1)=subplot(2,2,1)
edges=0:0.1:1;
[x y]=histc(osi_awake,edges)
%h1=hist(cell2mat(osi),edges)
b1=bar(edges,x,0.8,'facecolor', choose_color,'edgecolor',[0.28 0.28 0.28],'linewidth',0.5);
axis square
hold on
axis tight
ylim([0 4000])
plot(median(osi_awake), 0,'o','Color','r');
text(median(osi_awake),4000,sprintf('median = %.2f',median(osi_awake)),'Color',[0.8 0.2 0.2])
set(ax(1),'xminortick','on','tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(ax(1),'ticklength')*3)       
xlabel('OSI Awake','fontsize',14)
ylabel('# Boutons','fontsize',14)
        
% ----------------------
ax(2)=subplot(2,2,2)
edges=0:0.1:1;
[x y]=histc(dsi_awake,edges)
%h1=hist(cell2mat(osi),edges)
b1=bar(edges,x,0.8,'facecolor',choose_color,'edgecolor',[0.28 0.28 0.28],'linewidth',0.5);
axis square
hold on
plot(median(dsi_awake),0,'o','Color','r');
text(median(dsi_awake),4000,sprintf('median = %.2f',median(dsi_awake)),'Color',[0.8 0.2 0.2])
set(gca,'xminortick','on','tickdir','out')
axis tight
ylim([0 4000])
set(ax(2),'xminortick','on','tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(ax(2),'ticklength')*3)        
xlabel('DSI Awake','fontsize',14)
ylabel('# Boutons','fontsize',14)

ax(3)=subplot(2,2,3);
t = 0 : .01 : 2 * pi;
P = polar(t, 750 * ones(size(t)));
set(P, 'Visible', 'off');
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_ori_angle_awake;
r1=rose((2*pi/360)*ori_pref2,30);
set(r1,'color',choose_color);
set(gca,'xdir','reverse');

ax(4)=subplot(2,2,4)
t = 0 : .01 : 2 * pi;
P = polar(t, 500 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_dir_angle_awake;
r2=rose((2*pi/360)*ori_pref2,30);
set(r2,'color',choose_color);
set(gca,'xdir','reverse');

%%
set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'figure4_data_osi_dsi_37expt2.pdf']);
%% NORMALIZED
%%
dsi_awake= dsi_aw_norm;
osi_awake=osi_aw_norm;
pref_ori_angle_awake=prefosi_aw_norm(find(osi_awake>0.2));;
pref_dir_angle_awake=prefdsi_aw_norm(find(dsi_awake>0.2));
choose_color=[0 0 0];
%%
figure('name','norm awake')
clf
%awake_color=[110/255 150/255 200/255];
%anesthesia_color=[185/255 190/255 190/255];
awake_color=[0 0 0];
ax=[];
ax(1)=subplot(2,2,1)
edges=0:0.1:1;
[x y]=histc(osi_awake,edges)
%h1=hist(cell2mat(osi),edges)
b1=bar(edges,x,0.8,'facecolor', choose_color,'edgecolor',[0.28 0.28 0.28],'linewidth',0.5);
axis square
hold on
axis tight
ylim([0 4000])
plot(median(osi_awake), 0,'o','Color','r');
text(median(osi_awake),4000,sprintf('median = %.2f',median(osi_awake)),'Color',[0.8 0.2 0.2])
set(ax(1),'xminortick','on','tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(ax(1),'ticklength')*3)       
xlabel('OSI Awake','fontsize',14)
ylabel('# Boutons','fontsize',14)
        
% ----------------------
ax(2)=subplot(2,2,2)
edges=0:0.1:1;
[x y]=histc(dsi_awake,edges)
%h1=hist(cell2mat(osi),edges)
b1=bar(edges,x,0.8,'facecolor',choose_color,'edgecolor',[0.28 0.28 0.28],'linewidth',0.5);
axis square
hold on
plot(median(dsi_awake),0,'o','Color','r');
text(median(dsi_awake),4000,sprintf('median = %.2f',median(dsi_awake)),'Color',[0.8 0.2 0.2])
set(gca,'xminortick','on','tickdir','out')
axis tight
ylim([0 4000])
set(ax(2),'xminortick','on','tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',14,'ticklength',get(ax(2),'ticklength')*3)        
xlabel('DSI Awake','fontsize',14)
ylabel('# Boutons','fontsize',14)

ax(3)=subplot(2,2,3);
t = 0 : .01 : 2 * pi;
P = polar(t, 500 * ones(size(t)));
set(P, 'Visible', 'off');
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_ori_angle_awake;
r1=rose((2*pi/360)*ori_pref2,30);
set(r1,'color',choose_color);
set(gca,'xdir','reverse');

ax(4)=subplot(2,2,4)
t = 0 : .01 : 2 * pi;
P = polar(t, 200 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_dir_angle_awake;
r2=rose((2*pi/360)*ori_pref2,30);
set(r2,'color',choose_color);
set(gca,'xdir','reverse');

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'figure4_data_osi_dsi_37expt2_normalized.pdf']);
