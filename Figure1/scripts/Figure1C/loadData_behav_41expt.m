
filepathanalysis='G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour';
expt_finalName = [filepathanalysis, '\data_behav_41expt.mat'];
expt_finalName = [filepathanalysis, '\data_behav_40expt.mat'];
load(expt_finalName)
clear all
%% done on 03 July 2018
% plot behaviour for all data set all experiments 
newdir_behav={'160919_KS168_2P_KS\run01_behav_not_rand'
    '160916_KS166_2P_KS\run01_behav_not_rand'
    '160920_KS169_2P_KS\run01_behav_not_rand'
    '160920_KS168_2P_KS\run01_behav_not_rand'
    '160926_KS168_2P_KS\run01_behav_not_rand'
    '160926_KS169_2P_KS\run01_behav_not_rand'
    '160630_KS164_2P_KS\run03_ori12_V1_awake'
    '160630_KS167_2P_KS\run03_ori12_V1_awake'
    '160607_KS166_2P_KS\run03_ori12_V1'
    '160621_KS166_2P_KS\run03_ori12_V1_awake'
    '160703_KS164_2P_KS\run03_ori12_V1_awake'
    '160712_KS166_2P_KS\run03_ori12_V1_awake'
    '160712_KS167_2P_KS\run03_ori12_V1_awake'
    '160707_KS164_2P_KS\run03_ori12_V1_awake'
    '160913_KS169_2P_KS\run03_ori12_V1_awake'
    '160913_KS168_2P_KS\run03_ori12_V1_awake'
    %'160510_KS164_etl_2P_KS\run03_ori12_V1'% 100 pulses removed on 180705
    '160712_KS167_2P_KS\run03_ori12_V1_awake_post'
    '160712_KS167_2P_KS\run03_ori12_V1_awake_post2'
    '170110_KS174_2P_KS\run03_ori12_V1_awake'
    '170108_KS174_2P_KS\run03_ori12_V1_awake'
    '170110_KS173_2P_KS\run03_ori12_V1_awake'
    '160912_KS169_2P_KS\run03_ori12_V1_awake_post'
    '160913_KS169_2P_KS\run03_ori12_V1_awake_post'
    '140623_KS092_2P_KS\run01_ori_ds_V1'
    '140623_KS093_2P_KS\run01_ori_ds_V1'
    '150915_KS145_2P_KS\run01_ori_ds_V1_full'
    '140810_KS103_2P_KS\run01_ori_ds_V1_full'
    '151007_KS145_2P_KS\run012_ori_ds_V1_full'
    '150917_KS145_2P_KS\run02_ori_ds_V1_full'
    '140810_KS103_2P_KS\run02_ori_ds_V1_full'
    '151007_KS145_2P_KS\run02_ori_ds_V1_full'
    '151007_KS145_2P_KS\run023_ori_ds_V1_full'
    '140808_KS093_2P_KS\run03_ori_ds_V1_full'
    '150915_KS145_2P_KS\run034_ori_ds_V1_full'
    '140810_KS103_2P_KS\run03_ori_ds_V1_full'
    '151007_KS145_2P_KS\run034_ori_ds_V1_full'
    '151007_KS145_2P_KS\run03_ori_ds_V1_full'
    '170729_JC027_2P_JC\run02_ori12_V1'
    '170729_BV201_2P_JC\run01_ori12_L4'
    '170729_JC027_2P_JC\run01_ori12_V1_T6'};
%%
  %%
for ii=1:size(newdir_behav,1)
    size(newdir_behav,2)-1
    expt_behav{ii} = frGetExpt(newdir_behav{ii});
    
    expt2_behav{ii} =doLoadStimLogs3(expt_behav{ii});
    epochs_behav{ii} = expt2_behav{ii}.frames.epochs;
    stims_behav{ii}= expt2_behav{ii}.frames.stims;
    blanks_behav{ii}=expt2_behav{ii}.frames.blanks;
    stimulus_behav{ii}=[expt2_behav{ii}.prot.pars.ori];
    vel_behav{ii} = load([expt_behav{ii}.dirs.analrootpn,'\','velo_final.mat']);
    velocity_behav{ii} = vel_behav{ii}.velo_final(:,2);
    eye_behav{ii} = load([expt_behav{ii}.dirs.analrootpn,'\','diameter_data.mat']);
    
end
velocity_behav{16}(28443:end)=0;


%% save file data_behav_41expt.mat

% data_behav_40expt.diameter=data_behav_41expt.diameter(~cellfun('isempty',data_behav_41expt.diameter))  
% data_behav_40expt.velocity=data_behav_41expt.velocity(~cellfun('isempty',data_behav_41expt.velocity))  
% data_behav_40expt.newdir=data_behav_41expt.newdir(~cellfun('isempty',data_behav_41expt.newdir))  
% data_behav_40expt.stimulus=data_behav_41expt.stimulus(~cellfun('isempty',data_behav_41expt.stimulus))  
% data_behav_40expt.expt2=data_behav_41expt.expt2(~cellfun('isempty',data_behav_41expt.expt2))  
% data_behav_40expt.expt=data_behav_41expt.expt(~cellfun('isempty',data_behav_41expt.expt))  

data_behav_41expt.diameter=diameter_data_41expt;
data_behav_41expt.velocity=velocity_behav;
data_behav_41expt.newdir=newdir_behav;
data_behav_41expt.stimulus=stimulus_behav;
data_behav_41expt.expt2=expt2_behav;
data_behav_41expt.expt=expt_behav;
data_behav_41expt_v2=data_behav_41expt;

filepathanalysis='G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour';
expt_finalName = [filepathanalysis, '\data_behav_40expt.mat'];
%save(expt_finalName,'-v7.3','data_behav_40expt');

%% below not necessary - because this only loads data for analysis
% commented by KSocha 04-07-2023

% %%
% % figure (15)
% % clf
% % ii=40
% clear tr1
% clear tr2
% % % expt2=expt2_behav{ii}
% tr1=eye_behav{ii}.diameter_data.diameter_smooth_original;
% tr2=eye_behav{ii}.diameter_data.diameter_no_artifact;
% % % subplot(211)
% % % plot(tr1,'r');
% % % axis tight
% % % 
% % % subplot(212)
% % % plot(tr2,'k');
% % % axis tight
% % 
% % subplot(211);
% % plot(nanmean(reshape(tr1,[size(tr1,1)/expt2{ii}.info.nTrials expt2{ii}.info.nTrials]),2));
% % subplot(212)
% % plot(nanmean(reshape(tr2,[size(tr1,1)/expt2{ii}.info.nTrials expt2{ii}.info.nTrials]),2));
% %%
% diameter_data_41expt{ii}=tr1;
% diameter_data_41expt{ii}=tr2;
%% 

% data_behav_40expt.diameter=data_behav_41expt.diameter(~cellfun('isempty',data_behav_41expt.diameter))  
% data_behav_40expt.velocity=data_behav_41expt.velocity(~cellfun('isempty',data_behav_41expt.velocity))  
% data_behav_40expt.newdir=data_behav_41expt.newdir(~cellfun('isempty',data_behav_41expt.newdir))  
% data_behav_40expt.stimulus=data_behav_41expt.stimulus(~cellfun('isempty',data_behav_41expt.stimulus))  
% data_behav_40expt.expt2=data_behav_41expt.expt2(~cellfun('isempty',data_behav_41expt.expt2))  
% data_behav_40expt.expt=data_behav_41expt.expt(~cellfun('isempty',data_behav_41expt.expt))  

data_behav_41expt.diameter=diameter_data_41expt;
data_behav_41expt.velocity=velocity_behav;
data_behav_41expt.newdir=newdir_behav;
data_behav_41expt.stimulus=stimulus_behav;
data_behav_41expt.expt2=expt2_behav;
data_behav_41expt.expt=expt_behav;
data_behav_41expt_v2=data_behav_41expt;

filepathanalysis='G:\mousebox\code\mouselab\users\karolina\figure_paper_data_reanalyzed_behaviour';
expt_finalName = [filepathanalysis, '\data_behav_40expt.mat'];
%save(expt_finalName,'-v7.3','data_behav_40expt');

% subplot(211)
% plot(velocity_behav{ii});
% axis tight