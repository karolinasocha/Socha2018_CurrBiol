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
    filepathanalysis=[expt{ii}.dirs.analysis,'\',expt{ii}.animal,'\',expt{ii}.name];
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
nonoverlap_osi_awake=intersect(find(ori_awake>0.2),find(dir_awake<0.2))
nonoverlap_dsi_awake=intersect(find(ori_awake<0.2),find(dir_awake>0.2))

overlap_osi_dsi_awake=intersect(find(ori_awake>0.2),find(dir_awake>0.2))

%% select NON-OVERLAPPING

ori_pref_awake_nonoverlap=ori_pref_awake(nonoverlap_osi_awake);;
ori_pref_awake_overlap=ori_pref_awake(overlap_osi_dsi_awake);

dir_pref_awake_nonoverlap=dir_pref_awake(nonoverlap_dsi_awake);
dir_pref_awake_overlap=dir_pref_awake(overlap_osi_dsi_awake);

%% ANESTHESIA

%% PLOT FIGURE 4; qq plots of fraction of experiments
newdirs_anesthesia={'160707_KS164_2P_KS\run03_ori12_V1_anesthesia',...
    '160712_KS167_2P_KS\run03_ori12_V1_anesthesia',...
    '160621_KS166_2P_KS\run03_ori12_V1_anesthetized2',...
    '170110_KS173_2P_KS\run03_ori12_V1_anesthesia',...
    '170106_KS174_2P_KS\run03_ori12_V1_anesthesia',...
    '170110_KS174_2P_KS\run03_ori12_V1_anesthesia'}

clear stats_boutons_tmp
clear stats_boutons_normalized_tmp
for ii=1:size(newdirs_anesthesia,2)
 
    clear stats_version_mean_values
    expt{ii} = frGetExpt(newdirs_anesthesia{ii});
    expt2{ii} =doLoadStimLogs3(expt{ii});
    epochs{ii} = expt2{ii}.frames.epochs;
    stims{ii} = expt2{ii}.frames.stims;
    filepathanalysis=[expt{ii}.dirs.analysis,'\',expt{ii}.animal,'\',expt{ii}.name];
    fileNametcs=[filepathanalysis,'\stats_version_mean_values.mat'];
    load(fileNametcs)

stats_boutons_tmp{ii}=stats_version_mean_values.stats_boutons;
stats_boutons_normalized_tmp{ii}=stats_version_mean_values.stats_boutons_normalized;

end

clear stats_boutons_ansth
clear stats_boutons_normalized_ansth
for iAn=1:size(newdir,2)    
    stats_boutons_ansth{iAn}=stats_boutons_tmp{iAn};    
    stats_boutons_normalized_ansth{iAn}=stats_boutons_normalized_tmp{iAn};
end

%% transfer into arrays
stats_boutons_ansth=cell2mat(stats_boutons_ansth);
stats_boutons_normalized_ansth=cell2mat(stats_boutons_normalized_ansth);
%% orientation and direction tuning 
ori_ansth=[stats_boutons_ansth.ori_vector_sum_tune];
dir_ansth=[stats_boutons_ansth.dir_vector_sum_tune];
ori_norm_ansth=[stats_boutons_normalized_ansth.ori_vector_sum_tune];
dir_norm_ansth=[stats_boutons_normalized_ansth.dir_vector_sum_tune];

ori_pref_ansth=[stats_boutons_ansth.ori_vector_sum_angle_degree];
dir_pref_ansth=[stats_boutons_ansth.dir_vector_sum_angle_degree];

ori_pref_norm_ansth=[stats_boutons_normalized_ansth.ori_vector_sum_angle_degree];
dir_pref_norm_ansth=[stats_boutons_normalized_ansth.dir_vector_sum_angle_degree];
%% 
nonoverlap_osi_ansth=intersect(find(ori_ansth>0.2),find(dir_ansth<0.2))
nonoverlap_dsi_ansth=intersect(find(ori_ansth<0.2),find(dir_ansth>0.2))

overlap_osi_dsi_ansth=intersect(find(ori_ansth>0.2),find(dir_ansth>0.2))
%%
[h_dir,p_dir]=ranksum(dir_pref_ansth,dir_pref_awake) % p=1.5465e-09
[h_ori,p_ori]=ranksum(ori_pref_ansth,ori_pref_awake) % p=0.1428

%% select NON-OVERLAPPING


ori_pref_ansth_nonoverlap=ori_pref_ansth(nonoverlap_osi_ansth);
ori_pref_ansth_overlap=ori_pref_ansth(overlap_osi_dsi_ansth);

dir_pref_ansth_nonoverlap=dir_pref_ansth(nonoverlap_dsi_ansth);
dir_pref_ansth_overlap=dir_pref_ansth(overlap_osi_dsi_ansth);

%%
ori_pref_awake=ori_pref_awake_nonoverlap
ori_pref_anesthesia=ori_pref_ansth_nonoverlap
ranksum(ori_pref_awake,ori_pref_anesthesia )

ori_pref_awake=ori_pref_awake_overlap
ori_pref_anesthesia=ori_pref_ansth_overlap
ranksum(ori_pref_awake,ori_pref_anesthesia )

dir_pref_awake=dir_pref_awake_nonoverlap
dir_pref_anesthesia=dir_pref_ansth_nonoverlap
ranksum(dir_pref_awake,dir_pref_anesthesia)

dir_pref_awake=dir_pref_awake_overlap
dir_pref_anesthesia=dir_pref_ansth_overlap
ranksum(dir_pref_awake,dir_pref_anesthesia)
%%
[h_dir_overlap,p_dir_overlap]=ranksum(dir_pref_awake_overlap,dir_pref_ansth_overlap) % p=0.5846
[h_dir_nonoverlap,p_dir_nonoverlap]=ranksum(dir_pref_awake_nonoverlap,dir_pref_ansth_nonoverlap) % 0.1067
[h_ori_overlap,p_ori_overlap]=ranksum(ori_pref_awake_overlap,ori_pref_ansth_overlap) % p=0.3894
[h_ori_nonoverlap,p_ori_nonoverlap]=ranksum(ori_pref_awake_nonoverlap,ori_pref_ansth_nonoverlap) % p=0.7161
[h_ori,p_ori]=ranksum(ori_pref_ansth,ori_pref_awake) % p=0.1428
[h_dir,p_dir]=ranksum(dir_pref_ansth,dir_pref_awake) % p=1.5465e-09
%%
% Initialize arrays to store the p-values
p_dir_nonoverlaps = [];
p_ori_overlaps = [];
p_dir_overlaps = [];
p_ori_nonoverlaps = [];

% Loop through the iterations
for i = 1:1000
    numPoints =200;
    
    % Sample for dir_pref_awake_nonoverlap and dir_pref_ansth_nonoverlap
    sampled_awake = datasample(dir_pref_awake_nonoverlap, numPoints, 'Replace', true);
    sampled_ansth = datasample(dir_pref_ansth_nonoverlap, numPoints, 'Replace', true);

    [p_dir_nonoverlap, h_dir_nonoverlap] = ranksum(sampled_awake, sampled_ansth);
    % Append p_dir_nonoverlap to its array
    p_dir_nonoverlaps = [p_dir_nonoverlaps, p_dir_nonoverlap];
    
    % Sample for dir_pref_awake_overlap and dir_pref_ansth_overlap
    sampled_awake = datasample(dir_pref_awake_overlap, numPoints, 'Replace', true);
    sampled_ansth = datasample(dir_pref_ansth_overlap, numPoints, 'Replace', true);

    [p_dir_overlap, h_dir_overlap] = ranksum(sampled_awake, sampled_ansth);
    % Append p_dir_overlap to its array
    p_dir_overlaps = [p_dir_overlaps, p_dir_overlap];
    
    % Resample for dir_pref_awake_nonoverlap and dir_pref_ansth_nonoverlap
    sampled_awake = datasample(dir_pref_awake_nonoverlap, numPoints, 'Replace', true);
    sampled_ansth = datasample(dir_pref_ansth_nonoverlap, numPoints, 'Replace', true);

    [p_ori_nonoverlap, h_ori_nonoverlap] = ranksum(sampled_awake, sampled_ansth);
    % Append p_ori_nonoverlap to its array
    p_ori_nonoverlaps = [p_ori_nonoverlaps, p_ori_nonoverlap];
    
    % Sample for ori_pref_awake_overlap and ori_pref_ansth_overlap
    sampled_awake = datasample(ori_pref_awake_overlap, numPoints, 'Replace', true);
    sampled_ansth = datasample(ori_pref_ansth_overlap, numPoints, 'Replace', true);

    [p_ori_overlap, h_ori_overlap] = ranksum(sampled_awake, sampled_ansth);
    % Append p_ori_overlap to its array
    p_ori_overlaps = [p_ori_overlaps, p_ori_overlap];
end

prctile(p_ori_overlaps,50) % p=0.2431
prctile(p_ori_nonoverlaps,50) % p=0.0203
prctile(p_dir_overlaps,50) % p=0.4152
prctile(p_dir_nonoverlaps,50) % p=0.0198
%%
numPoints=1000
sampled_awake = datasample(ori_pref_awake_nonoverlap, numPoints, 'Replace', true);
sampled_ansth = datasample(ori_pref_ansth_nonoverlap, numPoints, 'Replace', true);
[p_ori_nonoverlap, h_ori_overlap] = ranksum(sampled_awake, sampled_ansth); %  0.3781

sampled_awake = datasample(ori_pref_awake_overlap, numPoints, 'Replace', true);
sampled_ansth = datasample(ori_pref_ansth_overlap, numPoints, 'Replace', true);
[p_ori_overlap, h_ori_overlap] = ranksum(sampled_awake, sampled_ansth); % 1.8059e-04


sampled_awake = datasample(dir_pref_awake_nonoverlap, numPoints, 'Replace', true);
sampled_ansth = datasample(dir_pref_ansth_nonoverlap, numPoints, 'Replace', true);
[p_dir_nonoverlap, h_dir_nonoverlap] = ranksum(sampled_awake, sampled_ansth); %  3.9707e-10

sampled_awake = datasample(dir_pref_awake_overlap, numPoints, 'Replace', true);
sampled_ansth = datasample(dir_pref_ansth_overlap, numPoints, 'Replace', true);
[p_dir_overlap, h_dir_overlap] = ranksum(sampled_awake, sampled_ansth); % 0.0963

%% COMPARISON WITH ENTIRE ANESTHESIA POPULATION
numPoints=1000
sampled_awake = datasample(ori_pref_awake_nonoverlap, numPoints, 'Replace', true);
sampled_ansth = datasample(ori_pref_ansth, numPoints, 'Replace', true);
[p_ori_nonoverlap, h_ori_nonoverlap] = ranksum(sampled_awake, sampled_ansth); %  0.3781

sampled_awake = datasample(ori_pref_awake_overlap, numPoints, 'Replace', true);
sampled_ansth = datasample(ori_pref_ansth, numPoints, 'Replace', true);
[p_ori_overlap, h_ori_overlap] = ranksum(sampled_awake, sampled_ansth); % 1.8059e-04

[pp,hh]= ranksum(ori_pref_awake_nonoverlap, ori_pref_ansth)% 0.1428
[pp,hh]= ranksum(ori_pref_awake_overlap, ori_pref_ansth) % 0.0862
[pp,hh]= ranksum(dir_pref_awake_overlap, dir_pref_ansth) % 0.0407
[pp,hh]= ranksum(dir_pref_awake_nonoverlap, dir_pref_ansth) %  1.5465e-09
%%
choose_color=[0 0 0];

%% get overlapping boutons:

pref_ori_angle_awake=ori_pref_awake_overlap;
pref_dir_angle_awake=dir_pref_awake_overlap;

figure('name','awake')
clf

ax(3)=subplot(2,2,3);
t = 0 : .01 : 2 * pi;
P = polar(t, 100 * ones(size(t)));
set(P, 'Visible', 'off');
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_ori_angle_awake;
r1=rose((2*pi/360)*ori_pref2,30);
set(r1,'color',choose_color);
set(gca,'xdir','reverse');

ax(4)=subplot(2,2,4)
t = 0 : .01 : 2 * pi;
P = polar(t, 50 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
%pref_osi_awake=pref_ori_angle_awake;
ori_pref2=pref_dir_angle_awake;
r2=rose((2*pi/360)*ori_pref2,30);
set(r2,'color',choose_color);
set(gca,'xdir','reverse');

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\Figure4A_overlapping_DS_OS_boutons\']; 
print(gcf,'-dpdf',[filepathanalysis, 'Figure4A_overlap_OS_DS_boutons_preference_awake_subset.pdf']);


%%
figure

ori_pref_awake=ori_pref_awake_overlap
ori_pref_anesthesia=ori_pref_ansth_overlap

ranksum(dir_pref_awake,dir_pref_anesthesia)
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

xlabel 'Awake Overlap';
ylabel 'Ansth Overlap';


% ------------------------------------
dir_pref_awake=dir_pref_awake_overlap
dir_pref_anesthesia=dir_pref_ansth_overlap

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

xlabel 'Awake Overlap';
ylabel 'Ansth Overlap';


ori_pref_awake=ori_pref_awake_nonoverlap
ori_pref_anesthesia=ori_pref_ansth_nonoverlap

ranksum(dir_pref_awake,dir_pref_anesthesia)
c=0;
a1 = angle(exp(j*(ori_pref_awake/180*pi+c)))/pi*180; % ANESTHESIA
a2 = angle(exp(j*(ori_pref_anesthesia/180*pi+c)))/pi*180; %AWAKE
%a1_norm = angle(exp(j*(ori_pref_norm_awake/180*pi+c)))/pi*180; % AWAKE 

v =0:10:100;
p1 = prctile(a1,v);
p2 = prctile(a2,v);
%p1_norm=prctile(a1_norm,v);

subplot(234)
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

xlabel 'Awake NoNoverlap';
ylabel 'Ansth NoNoverlap';


% ------------------------------------
dir_pref_awake=dir_pref_awake_nonoverlap
dir_pref_anesthesia=dir_pref_ansth_nonoverlap

c=pi/2
a1 = angle(exp(j*(dir_pref_awake/180*pi+c)))/pi*180; % AWAKE 
a2 = angle(exp(j*(dir_pref_anesthesia/180*pi+c)))/pi*180; % ANESTHESIA
%a1_norm = angle(exp(j*(dir_pref_norm_awake/180*pi+c)))/pi*180; % AWAKE 


v =0:10:100;
p1 = prctile(a1,v);
p2 = prctile(a2,v);
%p1_norm=prctile(a1_norm,v);
% plot(a1,dir_pref_awake,'.');
subplot(235)
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

xlabel 'Awake NoNoverlap';
ylabel 'Ansth NoNoverlap';

set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\Paper_figures2018\']; 
%print(gcf,'-dpdf',[filepathanalysis, 'figure4_data_qqplots_ansth_awake_norm_5expt.pdf']);
%%
[h_ansth,p_ansth]=kstest2(dir_pref_ansth_nonoverlap, dir_pref_ansth_overlap) % 0.0701
[h_awake,p_awake]=kstest2(dir_pref_awake_nonoverlap, dir_pref_awake_overlap)

[h_ansth_ori,p_ansth_ori]=kstest2(ori_pref_ansth_nonoverlap, ori_pref_ansth_overlap) % 0.0701
[h_awake_ori,p_awake_ori]=kstest2(ori_pref_awake_nonoverlap, ori_pref_awake_overlap)

[h_awake_ansth_non_ori,p_awake_ansth_non_ori]=kstest2(ori_pref_awake_nonoverlap, ori_pref_ansth_nonoverlap)
[h_awake_ansth_over_ori,p_awake_ansth_over_ori]=kstest2(ori_pref_awake_overlap, ori_pref_ansth_overlap)

[h_awake_ansth_non_dir,p_awake_ansth_non_dir]=kstest2(dir_pref_awake_nonoverlap, dir_pref_ansth_nonoverlap)
[h_awake_ansth_over_dir,p_awake_ansth_over_dir]=kstest2(dir_pref_awake_overlap, dir_pref_ansth_overlap)




