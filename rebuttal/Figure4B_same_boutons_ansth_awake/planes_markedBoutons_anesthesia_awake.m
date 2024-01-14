clear all

newdirs_anesthesia={'160707_KS164_2P_KS\run03_ori12_V1_anesthesia',...
    '160712_KS167_2P_KS\run03_ori12_V1_anesthesia',...
    '160621_KS166_2P_KS\run03_ori12_V1_anesthetized2',...
    '170110_KS173_2P_KS\run03_ori12_V1_anesthesia',...
    '170106_KS174_2P_KS\run03_ori12_V1_anesthesia',...
    '170110_KS174_2P_KS\run03_ori12_V1_anesthesia'}

newdirs_awake={'160707_KS164_2P_KS\run03_ori12_V1_awake',...
    '160712_KS167_2P_KS\run03_ori12_V1_awake',...
    '160621_KS166_2P_KS\run03_ori12_V1_awake',...
    '170110_KS173_2P_KS\run03_ori12_V1_awake',...
    '170106_KS174_2P_KS\run03_ori12_V1',...
    '170110_KS174_2P_KS\run03_ori12_V1_awake'};

%%
filepath='G:\mousebox\analysis\';

for i=2:6

newdir_anesthesia=newdirs_anesthesia{i};
newdir_awake=newdirs_awake{i};

clear tcs_handseg_ansth
clear tcs_handseg_awake

load(fullfile(filepath,newdir_anesthesia,'tcs_handseg.mat'));
% load('G:\mousebox\analysis\160707_KS164_2P_KS\run03_ori12_V1_anesthesia\tcs_handseg.mat')
tcs_handseg_ansth=tcs_handseg;
clear tcs_handseg
% load('G:\mousebox\analysis\160707_KS164_2P_KS\run03_ori12_V1_awake\tcs_handseg.mat')
load(fullfile(filepath,newdir_awake,'tcs_handseg.mat'))
tcs_handseg_awake=tcs_handseg;
%

if isfield(tcs_handseg_ansth, 'colormap') & isfield(tcs_handseg_awake, 'colormap')
    
    % Assuming tcs_handseg.colormap and tcs_handseg.mask are your images
    figure;
    clf
    subplot('Position',[0.05, 0.2, 0.45, 0.45])
    % Display the first image
    % imshow(tcs_handseg_ansth.colormap);
    imshow(tcs_handseg_ansth.colormap(:,:,3));
    % Hold on to overlay the second image
    hold on;
    % Display the second image with transparency
    h = imshow(tcs_handseg_ansth.mask);
    set(h, 'AlphaData', 0.3); % Adjust the transparency level as needed
    % Release the hold to prevent further overlay
    hold off;
    axis square 
    title('ANESTHESIA')

    subplot('Position',[0.55, 0.2, 0.45, 0.45])
    % Display the first image
    % imagesc(tcs_handseg_awake.colormap(:,:,1));
    imshow(tcs_handseg_awake.colormap(:,:,3));
    % Hold on to overlay the second image
    hold on;
    % Display the second image with transparency
    h = imshow(tcs_handseg_awake.mask);

    set(h, 'AlphaData', 0.3); % Adjust the transparency level as needed
    % Release the hold to prevent further overlay
    hold off;
    axis square 
    title('AWAKE')

    set(gcf,'paperunits','centimeters','papersize' ,[20,20],'color','w','paperposition',[0,0,20,20],'inverthardcopy','off');
    filepathanalysis=['G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\Figure4B_same_boutons_ansth_awake\']; 
    filename=sprintf(['Figure4B_samePlane_anesthesia_Awake_',newdir_awake(1:18),'.pdf'])

    print(gcf,'-dpdf',[filepathanalysis, filename]);

else
    disp('Field ''colormap'' does not exist in tcs_handseg_awake.');
end
end
%%




