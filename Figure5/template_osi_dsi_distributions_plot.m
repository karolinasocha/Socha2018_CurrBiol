%% PLOT FIGURES FOR SUPPLEMENTARY FIGURE 5
% corr plot preferences 
% load Data corr preferences
% 'G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences'
%% load first script with 23 experiments:  loadData_corr_preferences_sameanimals
%% PLOT SCATTER PLOTS FOR PAPER
figure(6)
clf
subplot(231);
clear colors_dots
clear colors_stimulus
clear color_boutons
scatter(ori_index_example,corr_bouton_pupil,2,[0 0 0],'filled');
hold on
plot([0.2 0.2],[-0.6 1],'--','color',[0.5 0.5 0.5]);
%
axis tight
xlabel('OSI','fontsize',12)
ylabel('Corr Bouton-Pupil','fontsize',12)
axis square
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],...
'xlim',[0 1],'xtick',[0:0.2:1],...
    'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4)       
title('Example experiment')


subplot(232);
clear colors_dots
clear colors_stimulus
clear color_boutons
scatter(ori_pref_example,corr_bouton_pupil,2,[0 0 0],'filled');
hold on
%list=find(dir_index_example>0.2);
%scatter(dir_pref_example(list),corr_bouton_pupil(list),5,[ 1 0 0],'filled');
hold on
axis tight
xlabel('OSI Pref','fontsize',12)
ylabel('Corr Bouton-Pupil','fontsize',12)
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 180],'xtick',[0:30:180],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4)       
title('Example experiment')
axis square
hold on
clear m_val_dir

ed = 0:22.5:180;
[ed, y_mu, y_s] = binSamples(ori_pref_example, corr_bouton_pupil, ed);
m_val_dir(iAn,:)=y_mu;
ed_dir=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5],'linewidth',2);


clear m_val_dir
for iAn=1:23
    
dir_pref_all=[stats_boutons{iAn}.dir_vector_sum_angle_degree];
ori_pref_all=[stats_boutons{iAn}.ori_vector_sum_angle_degree];
ori_index_all=[stats_boutons{iAn}.ori_vector_sum_tune];
dir_index_all=[stats_boutons{iAn}.dir_vector_sum_tune];
corr_bouton_pupil_all=correlationcoef{iAn}.corr_bouton_pupil;
corr_bouton_population_all=correlationcoef{iAn}.corr_bouton_population;

subplot(233);
clear ed
ed = 0:22.5:180;
clear y_mu
%list=find(ori_index_all>tresh_sel)
[ed, y_mu, y_s] = binSamples(ori_pref_all, corr_bouton_pupil_all, ed);
m_val_ori_pupil(iAn,:)=y_mu;
ed_dir=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(234);
clear y_mu
clear ed
ed = 0:0.1:1;
[ed, y_mu, y_s] = binSamples(ori_index_all, corr_bouton_pupil_all, ed);
m_val_ori_index_pupil(iAn,:)=y_mu;
ed_ori_index=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(235);
clear y_mu
clear ed
ed = 0:0.1:1;
[ed, y_mu, y_s] = binSamples(ori_index_all, corr_bouton_population_all, ed);
m_val_ori_index_pop(iAn,:)=y_mu;
ed_ori_index=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(236);
clear y_mu
clear ed
ed = 0:22.5:180;
[ed, y_mu, y_s] = binSamples(ori_pref_all, corr_bouton_population_all, ed);
m_val_ori_pop(iAn,:)=y_mu;
ed_ori=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

end

subplot(233);
hold on
plot(ed_dir,nanmean(m_val_ori_pupil,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 180],'xtick',[0:30:180],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('OSI Pref','fontsize',12);
title('All sessions')
ylabel('Corr','fontsize',12);
axis square

subplot(234);
hold on
plot(ed_ori_index,nanmean(m_val_ori_index_pupil,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 1],'xtick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('OSI','fontsize',12);
title('All sessions')
ylabel('Corr Pupil','fontsize',12);
axis square

subplot(235);
hold on
plot(ed_ori_index,nanmean(m_val_ori_index_pop,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 1],'xtick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('OSI','fontsize',12);
title('All sessions')
ylabel('Corr Population','fontsize',12);
axis square

subplot(236);
hold on
plot(ed_ori,nanmean(m_val_ori_pop,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 180],'xtick',[0:30:180],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('OSI Pref','fontsize',12);
title('All sessions')
ylabel('Corr Population','fontsize',12);
axis square

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences\']; 
print(gcf,'-dpdf',[filepathanalysis, '\fig5_SuppMat_corr_example_allsesions_OSI_paper.pdf']);

%% DIRECTION SM SUMMARY
figure(7)
clf

clear m_val_dir
for iAn=1:23
    
dir_pref_all=[stats_boutons{iAn}.dir_vector_sum_angle_degree];
ori_pref_all=[stats_boutons{iAn}.ori_vector_sum_angle_degree];
ori_index_all=[stats_boutons{iAn}.ori_vector_sum_tune];
dir_index_all=[stats_boutons{iAn}.dir_vector_sum_tune];
corr_bouton_pupil_all=correlationcoef{iAn}.corr_bouton_pupil;
corr_bouton_population_all=correlationcoef{iAn}.corr_bouton_population;

subplot(234);
clear ed
ed = 0:22.5:360;
clear y_mu
%list=find(ori_index_all>tresh_sel)
[ed, y_mu, y_s] = binSamples(dir_pref_all, corr_bouton_pupil_all, ed);
m_val_dir_pupil(iAn,:)=y_mu;
ed_dir=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(232);
clear y_mu
clear ed
ed = 0:0.1:1;
[ed, y_mu, y_s] = binSamples(dir_index_all, corr_bouton_pupil_all, ed);
m_val_dir_index_pupil(iAn,:)=y_mu;
ed_dir_index=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(233);
clear y_mu
clear ed
ed = 0:0.1:1;
[ed, y_mu, y_s] = binSamples(dir_index_all, corr_bouton_population_all, ed);
m_val_dir_index_pop(iAn,:)=y_mu;
ed_dir_index=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

subplot(231);
clear y_mu
clear ed
ed = 0:22.5:360;
[ed, y_mu, y_s] = binSamples(dir_pref_all, corr_bouton_population_all, ed);
m_val_dir_pop(iAn,:)=y_mu;
ed_dir=ed;
hold on
plot(ed,y_mu,'color',[0.5 0.5 0.5]);

end

subplot(234);
hold on
plot(ed_dir,nanmean(m_val_dir_pupil,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 360],'xtick',[0:30:360],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('DSI Pref','fontsize',12);
title('All sessions')
ylabel('Corr Pupil','fontsize',12);
axis square

subplot(232);
hold on
plot(ed_dir_index,nanmean(m_val_dir_index_pupil,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 1],'xtick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('DSI','fontsize',12);
title('All sessions')
ylabel('Corr Pupil','fontsize',12);
axis square

subplot(233);
hold on
plot(ed_dir_index,nanmean(m_val_dir_index_pop,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 1],'xtick',[0:0.2:1],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('DSI','fontsize',12);
title('All sessions')
ylabel('Corr Population','fontsize',12);
axis square

subplot(231);
hold on
plot(ed_dir,nanmean(m_val_dir_pop,1),'color','r','linewidth',2);

axis tight;
set(gca,'ylim',[-0.6 1],'ytick',[-0.6:0.2:1],'xlim',[0 360],'xtick',[0:30:360],'tickdir','out','box','off','tickdir','out',...
'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*4);      
xlabel('DSI Pref','fontsize',12);
title('All sessions')
ylabel('Corr Population','fontsize',12);
axis square

set(gcf,'paperunits','centimeters','papersize' ,[30,30],'color','w','paperposition',[0,0,30,30],'inverthardcopy','off')
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\figure_relation_corr_preferences\']; 
%print(gcf,'-dpdf',[filepathanalysis, '\fig5_SuppMat_corr_example_allsesions_DSI_paper.pdf']);


