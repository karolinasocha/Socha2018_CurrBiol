%% pupil size with differen sizes plots

clear all

%%
 newdir={'190217_XH315_behav_JC\run00_direction_cardinal_4Hz_sizes',...
 '190217_XH158_behav_JC\run00_direction_cardinal_4Hz_sizes',...
 '190217_XH322_behav_JC\run00_direction_cardinal_4Hz_sizes'}

% newdir={'190217_XH315_behav_JC\run00_direction_cardinal_4Hz_sizes',...
% '190217_XH158_behav_JC\run00_direction_cardinal_4Hz_sizes',...
% '190217_XH309_behav_JC\run00_direction_cardinal_4Hz_sizes',...
% '190217_XH322_behav_JC\run00_direction_cardinal_4Hz_sizes'}

for iAn=1:size(newdir,2)
clear tt
clear pupildata
filepathanalysis=['J:\data\analysis\' newdir{iAn}]; 
pupil_name = [filepathanalysis, '\pupildata.mat'];
load(pupil_name);
diameter_pupil{iAn}=pupildata.diameter;
diameter_stim{iAn}=pupildata.diameter.diam_stim;
diameter_epochs{iAn}=pupildata.diameter.diam_epochs
baseline_median{iAn}=nanmedian(pupildata.diameter.diameter);
diameter_av_epochs{iAn}=squeeze(nanmean(diameter_epochs{iAn},2));

clear X
    X=pupildata.diameter.diameter;
for istims=1:size(diameter_av_epochs{iAn},2)
    clear tmp2
    tmp2=diameter_av_epochs{iAn}(:,istims); % this is previous calculation in mm
    diameter_epochs_zscored{iAn}(:,istims) = bsxfun(@rdivide, bsxfun(@minus, tmp2, nanmean(X)) , ...
                     nanstd(X));
end
end

 %% interpolate pupil   

clear diam_interp
for iAn=1:size(newdir,2)
    for iStim=1:24
    clear dtime
    clear p2time
    dtime = linspace(0,1,size(diameter_av_epochs{iAn}(:,iStim),1));
    p2time =  linspace(0,1,315);
    diam_interp{iAn}(:,iStim) = interp1(dtime,diameter_av_epochs{iAn}(:,iStim),p2time);
    diam_interp_zscored{iAn}(:,iStim) = interp1(dtime,diameter_epochs_zscored{iAn}(:,iStim),p2time);
    
    end
end

%%
diameter_input=diam_interp_zscored;

%%
clear diameter_mean
stim0deg=[1 5 9 13 17 21];
stim90deg=1+stim0deg;
stim180deg=2+stim0deg;
stim270deg=3+stim0deg;

for iAn=1:size(newdir,2)
diam_stimulus0deg(iAn,:,:)=diameter_input{iAn}(:,stim0deg);
diam_stimulus90deg(iAn,:,:)=diameter_input{iAn}(:,stim90deg);
diam_stimulus180deg(iAn,:,:)=diameter_input{iAn}(:,stim180deg);
diam_stimulus270deg(iAn,:,:)=diameter_input{iAn}(:,stim270deg);
end

av_diam_stimulus0deg=squeeze(nanmean(diam_stimulus0deg,1));
av_diam_stimulus90deg=squeeze(nanmean(diam_stimulus90deg,1));
av_diam_stimulus180deg=squeeze(nanmean(diam_stimulus180deg,1));
av_diam_stimulus270deg=squeeze(nanmean(diam_stimulus270deg,1));

%%
av_diam_stimulus0deg=squeeze(nanmean(diam_stimulus0deg,1));
av_diam_stimulus90deg=squeeze(nanmean(diam_stimulus90deg,3));
av_diam_stimulus180deg=squeeze(nanmean(diam_stimulus180deg,3));
av_diam_stimulus270deg=squeeze(nanmean(diam_stimulus270deg,3));
%% calculate mean across mice +- SEM
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
%
title_str={'Full','90deg','70deg','50deg','30deg','10deg'}
clf
for iStim=1:size(diam_stimulus0deg,3)
ax_pupil_response(iStim) = axes('position',[0.15*(iStim-1)+.1,.45,.1,.2],'units','normalized');
end

colors_line=parula(12);

for iStim=1:size(diam_stimulus0deg,3);
mean_pupil_delta_awake=squeeze(diam_stimulus0deg(:,:,iStim));
ylimitation_pupil=[-1 1];

% Calculate the mean and standard deviation across columns
meanData_awake = nanmean(mean_pupil_delta_awake,1);
stdData_awake = nanstd(mean_pupil_delta_awake,1);
x = 1:numel(meanData_awake);

% Calculate the Standard Error of the Mean (SEM)
semData_awake = stdData_awake ./ sqrt(size(mean_pupil_delta_awake, 1));

axes(ax_pupil_response(iStim))
plot(x, meanData_awake, 'color',colors_line(iStim,:),'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_awake - semData_awake, fliplr(meanData_awake + semData_awake)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil 0 deg stim (zscored)', 'FontSize',16,'Color','k');
xlabel('Time (s)','FontSize',16,'Color','k');
%set(gca,'ylim',[-1, 1], 'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gca,'xtick',[63 221],'xticklabel',{'0' '5'},...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(gca,'ticklength')*8);

hold on
plot([63 63],ylimitation_pupil,'--k');
hold on
plot([221 221],ylimitation_pupil,'--k');
title(title_str{iStim})
end
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\SuppFigure1\'];

print(gcf,'-dpdf',[filepathanalysis, 'SupplFigure1B_3mice_differentSize_0deg_SEM_animals.pdf']);
%
%% 180 deg

%% calculate mean across mice +- SEM
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
%
title_str={'Full','90deg','70deg','50deg','30deg','10deg'}
clf
for iStim=1:size(diam_stimulus0deg,3)
ax_pupil_response(iStim) = axes('position',[0.15*(iStim-1)+.1,.45,.1,.2],'units','normalized');
end

colors_line=parula(12);

for iStim=1:size(diam_stimulus0deg,3);
mean_pupil_delta_awake=squeeze(diam_stimulus180deg(:,:,iStim));
ylimitation_pupil=[-1 1];

% Calculate the mean and standard deviation across columns
meanData_awake = nanmean(mean_pupil_delta_awake,1);
stdData_awake = nanstd(mean_pupil_delta_awake,1);
x = 1:numel(meanData_awake);

% Calculate the Standard Error of the Mean (SEM)
semData_awake = stdData_awake ./ sqrt(size(mean_pupil_delta_awake, 1));

axes(ax_pupil_response(iStim))
plot(x, meanData_awake, 'color',colors_line(iStim,:),'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_awake - semData_awake, fliplr(meanData_awake + semData_awake)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil 180 deg stim (zscored)', 'FontSize',16,'Color','k');
xlabel('Time (s)','FontSize',16,'Color','k');
%set(gca,'ylim',[-1, 1], 'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gca,'xtick',[63 221],'xticklabel',{'0' '5'},...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(gca,'ticklength')*8);

hold on
plot([63 63],ylimitation_pupil,'--k');
hold on
plot([221 221],ylimitation_pupil,'--k');
title(title_str{iStim})
end
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\SuppFigure1\'];

print(gcf,'-dpdf',[filepathanalysis, 'SupplFigure1B_3mice_differentSize_180deg_SEM_animals.pdf']);
%
%% 90 deg

%% calculate mean across mice +- SEM
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
%
title_str={'Full','90deg','70deg','50deg','30deg','10deg'}
clf
for iStim=1:size(diam_stimulus90deg,3)
ax_pupil_response(iStim) = axes('position',[0.15*(iStim-1)+.1,.45,.1,.2],'units','normalized');
end

colors_line=parula(12);

for iStim=1:size(diam_stimulus90deg,3);
mean_pupil_delta_awake=squeeze(diam_stimulus90deg(:,:,iStim));
ylimitation_pupil=[-1 1];

% Calculate the mean and standard deviation across columns
meanData_awake = nanmean(mean_pupil_delta_awake,1);
stdData_awake = nanstd(mean_pupil_delta_awake,1);
x = 1:numel(meanData_awake);

% Calculate the Standard Error of the Mean (SEM)
semData_awake = stdData_awake ./ sqrt(size(mean_pupil_delta_awake, 1));

axes(ax_pupil_response(iStim))
plot(x, meanData_awake, 'color',colors_line(iStim,:),'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_awake - semData_awake, fliplr(meanData_awake + semData_awake)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil 90 deg stim (zscored)', 'FontSize',16,'Color','k');
xlabel('Time (s)','FontSize',16,'Color','k');
%set(gca,'ylim',[-1, 1], 'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gca,'xtick',[63 221],'xticklabel',{'0' '5'},...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(gca,'ticklength')*8);

hold on
plot([63 63],ylimitation_pupil,'--k');
hold on
plot([221 221],ylimitation_pupil,'--k');
title(title_str{iStim})
end
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\SuppFigure1\'];

print(gcf,'-dpdf',[filepathanalysis, 'SupplFigure1B_3mice_differentSize_90deg_SEM_animals.pdf']);
%
%%
%% calculate mean across mice +- SEM
fig = figure('name',sprintf('Pupil size behaviour'),'color','w','paperunits',...
    'centimeters','papersize',[21,29.7],'paperposition',[0,0,21,29.7]);
%
title_str={'Full','90deg','70deg','50deg','30deg','10deg'}
clf
for iStim=1:size(diam_stimulus270deg,3)
ax_pupil_response(iStim) = axes('position',[0.15*(iStim-1)+.1,.45,.1,.2],'units','normalized');
end

colors_line=parula(12);

for iStim=1:size(diam_stimulus270deg,3);
mean_pupil_delta_awake=squeeze(diam_stimulus270deg(:,:,iStim));
ylimitation_pupil=[-1 1];

% Calculate the mean and standard deviation across columns
meanData_awake = nanmean(mean_pupil_delta_awake,1);
stdData_awake = nanstd(mean_pupil_delta_awake,1);
x = 1:numel(meanData_awake);

% Calculate the Standard Error of the Mean (SEM)
semData_awake = stdData_awake ./ sqrt(size(mean_pupil_delta_awake, 1));

axes(ax_pupil_response(iStim))
plot(x, meanData_awake, 'color',colors_line(iStim,:),'LineWidth', 2);
hold on;

% Create a shaded region for the SEM
fill([x, fliplr(x)], [meanData_awake - semData_awake, fliplr(meanData_awake + semData_awake)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
axis tight

ylabel('Pupil 270 deg stim (zscored)', 'FontSize',16,'Color','k');
xlabel('Time (s)','FontSize',16,'Color','k');
%set(gca,'ylim',[-1, 1], 'tickdir','out','fontsize',14,'ticklength',get(gca,'ticklength')*4);
box off

set(gca,'xtick',[63 221],'xticklabel',{'0' '5'},...
    'tickdir','out','box','off','layer','top','color','none',...
    'fontsize',14,'ticklength',get(gca,'ticklength')*8);

hold on
plot([63 63],ylimitation_pupil,'--k');
hold on
plot([221 221],ylimitation_pupil,'--k');
title(title_str{iStim})
end
set(gcf,'paperunits','centimeters','papersize' ,[21,29.7],'color','w','paperposition',[0,0,21,29.7],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\SuppFigure1\'];

print(gcf,'-dpdf',[filepathanalysis, 'SupplFigure1B_3mice_differentSize_270deg_SEM_animals.pdf']);
%