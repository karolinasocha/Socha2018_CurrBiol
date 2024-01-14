clear all
% this load processed data for checking statistical tests

% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_0_05.mat');
% meanData_0_05=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_05_1.mat');
meanData_05_1=meanData

load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_1_15.mat');
meanData_1_15=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_15_2.mat');
meanData_15_2=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_2_25.mat');
meanData_2_25=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_25_3.mat');
meanData_25_3=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_3_35.mat');
meanData_3_35=meanData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_35_4.mat');
meanData_35_4=meanData

load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\meanData_45_5.mat');
meanData_45_5=meanData
mean_Data=[meanData_05_1;meanData_1_15;meanData_15_2;meanData_2_25;meanData_25_3;meanData_3_35;meanData_35_4;;meanData_45_5];


%%
% load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_0_05.mat');
% semData_0_05=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_05_1.mat');
semData_05_1=semData

load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_1_15.mat');
semData_1_15=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_15_2.mat');
semData_15_2=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_2_25.mat');
semData_2_25=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_25_3.mat');
semData_25_3=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_3_35.mat');
semData_3_35=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_35_4.mat');
semData_35_4=semData
load('G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\semData_45_5.mat');
semData_45_5=semData

sem_Data=[semData_05_1;semData_1_15;semData_15_2;semData_2_25;semData_25_3;semData_3_35;semData_35_4;semData_45_5];


%%
time_calc=[1,1.5,2,2.5,3,3.5,4,4.5];
figure
for iDir=1:12
subplot('Position',[iDir*0.07+0.05,.4,0.05,0.3]);
% er=errorbar(time_calc, mean_Data(:,iDir), sem_Data(:,iDir), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% er.Color=[0,0,0];
% set(gca,'xtick',[2,3,4,4.5],'xticklabel',[2,3,4,5],'ylim',[-10, 30],'tickdir','out','box','off','tickdir','out',...
% 'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);
%    
er=errorbar(time_calc, mean_Data(:,iDir), sem_Data(:,iDir), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
er.Color=[0,0,0];

if iDir > 1
    % For subplots other than the first one, remove x and y tick labels
    set(gca,'xtick',[2,3,4,4.5],'xticklabel',[2,3,4,5],'ylim',[-10, 30],'tickdir','out','box','off','tickdir','out',...
    'layer','top','color','none','fontsize',10,'ticklength',get(gca,'ticklength')*3);
    set(gca, 'xticklabel', [], 'yticklabel', []);

else
    % For the first subplot, set x and y tick labels
    set(gca, 'xtick', [1, 2, 3, 4, 4.5], 'xticklabel', [1, 2, 3, 4, 5], ...
        'ylim', [-10, 30], 'tickdir', 'out', 'box', 'off', 'tickdir', 'out', ...
        'layer', 'top', 'color', 'none', 'fontsize', 10, 'ticklength', get(gca, 'ticklength') * 3);
    ylabel('delta pupil (%)','fontsize',14);
    xlabel('time from onset','fontsize',14);
end
end
% ax=axes('Position',[1,2,1,3]);
set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\']; 
print(gcf,'-dpdf',[filepathanalysis, 'Figure1E_pupil_effect_timeStimulation.pdf']);

%%
time_calc=[1.5,2,2.5,3,3.5,4,4.5];
figure;
for iDir=1:12
% Plot data points with error bars
subplot(3,4,iDir)
errorbar(time_calc, mean_Data(:,iDir), sem_Data(:,iDir), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
set(gca,'xticklabel',[1.5,2,2.5,3,3.5,4,5],'ylim',[-10, 30])

end
%%
iDir=1
x=time_calc;
y=mean_Data(:,iDir)';
coefficients = polyfit(x, y, 1);

% Generate points along the fitted line for plotting
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);

% Plot the original data
figure;
plot(x, y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;

% Plot the linear regression line
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);



