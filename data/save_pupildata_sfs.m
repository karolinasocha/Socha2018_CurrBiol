% newdirs={'240314_AE01__no_cam_SH', '240305_MG04__no_cam_SH','240311_SH056__no_cam_SH'}
% sessions={'run00_140423_12dir-sf08-tf4SinSeq5s-20T, 'run00_195835_12dir-sf08-tf-4SinSeq''}
% Experiments list

% EFFECT IN FEMALE MICE
% '240314_AE01__no_cam_SH\run00_140423_12dir-sf08-tf4SinSeq5s-20T' % % -FEMALE; sequence; 10 trials; poor quality to remove
% '240305_MG04__no_cam_SH\run00_195835_12dir-sf08-tf-4SinSeq' % -FEMALE; % sequence; 10 trials; no effect

% EFFECT DIFFERENT DURATION
% '240311_SH056__no_cam_SH\run00_154140_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH056__no_cam_SH\run01_160616_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH056__no_cam_SH\run02_162601_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH056__no_cam_SH\run03_164335_12dir-sf08-tf4SinSeq2s' - 2s duration

% '240311_SH057__no_cam_SH\run00_170314_12dir-sf08-tf4SinSeq5s' - 5s duration
% '240311_SH057__no_cam_SH\run01_172413_12dir-sf08-tf4SinSeq4s' - 4s duration
% '240311_SH057__no_cam_SH\run02_174356_12dir-sf08-tf4SinSeq3s' - 3s duration
% '240311_SH057__no_cam_SH\run03_180209_12dir-sf08-tf4SinSeq2s' - 2s duration

% EFFECT DIFFERENT TFs
% '240313_SH055__no_cam_SH\run00_133318_5tf-4dir-sf08-SinRand-5s' 
% '240313_SH056__no_cam_SH\run00_145132_5tf-4dir-sf08-SinRand-5s'
% '240313_SH057__no_cam_SH\run00_160812_5tf-4dir-sf08-SinRand-5s'

% EFFECT DIFFERENT SFs
% '240313_SH055__no_cam_SH\run01_140753_5sf-4dir-tf4-SinRand5s'
% '240313_SH056__no_cam_SH\run01_152943_5sf-4dir-tf4-SinRand5s'
% '240313_SH057__no_cam_SH\run01_164309_5sf-4dir-tf4-SinRand5s'

%% read mat file with pupil trace and timestamps

clear all

path_data='G:\mousebox\code\mouselab\users\karolina\FiguresPaper2023\rebuttal\data_rebuttal\data2\'
newdir='240313_SH055__no_cam_SH'
session='run01_140753_5sf-4dir-tf4-SinRand5s'
filepathanalysis=fullfile(path_data,  newdir, session); 
expt_name = [filepathanalysis, '\KSdata.mat'];
% 'SH055_12dir_SF0.08_TF4_5sstim_5sblank_Seq_20tr'
experiment_data=load(expt_name);
fps=experiment_data.fsl
stimuli_params=experiment_data.stimuli
pupil_data=experiment_data.timepar
pupil_timestamps=pupil_data(:,1);
pupil_area=pupil_data(:,2);
pupil_area_zscore = bsxfun(@rdivide, bsxfun(@minus, pupil_area, nanmean(pupil_area)) , nanstd(pupil_area));
pupil_area_zscore2=zscore(pupil_area);
stimulusOnset=experiment_data.vstimon;

%%

%%
stim_directions=unique(stimuli_params(:,1));
stim_sfs=unique(stimuli_params(:,2));
ntrials= length(intersect(find(stimuli_params(:,1)==0),find(stimuli_params(:,2)==0.08 )));

stim_duration=5;
prestim_duration=5;
poststim_duration=5;

pupil_trials_stim=nan(length(stim_directions),length(stim_sfs), ntrials,22*stim_duration); % 110 for 5s
pupil_trials_blanks=nan(length(stim_directions),length(stim_sfs), ntrials,22*stim_duration); % 110 for 5s
pupil_trials_poststim=nan(length(stim_directions),length(stim_sfs), ntrials,22*stim_duration);
pupil_trials_epochs=nan(length(stim_directions),length(stim_sfs), ntrials,110+22*stim_duration); % 220 for 5s
pupil_trials_entire_epochs=nan(length(stim_directions),length(stim_sfs), ntrials,220+22*stim_duration); % 330 for 5s

for istim=1:length(stim_directions)
    for isf = 1:length(stim_sfs)
        
        indx=intersect(find(stimuli_params(:,1)==stim_directions(istim)),find(stimuli_params(:,2)==stim_sfs(isf)));
        stimOnset_indx=stimulusOnset(indx);
        stimOffset_indx=round(stim_duration*fps)+stimOnset_indx; 
        stimOffset_indx_fixed=round(5*fps)+stimOnset_indx; 

            for itrial=1:length(stimOnset_indx)

                target_timestamp = stimOnset_indx(itrial);
                differences = abs(pupil_timestamps - target_timestamp);
                [~, nearestIndex_onset] = min(differences);

                target_timestampOffset = stimOffset_indx(itrial);
                differences2 = abs(pupil_timestamps - target_timestampOffset);
                [~, nearestIndex_offset] = min(differences2); 

                target_timestamppre = stimOnset_indx(itrial)-prestim_duration*fps;
                differences3 = abs(pupil_timestamps - target_timestamppre);
                [~, nearestIndex_onsetPre] = min(differences3);

                target_timestamppost = stimOffset_indx(itrial)+poststim_duration*fps;
                differences4 = abs(pupil_timestamps - target_timestamppost);
                [~, nearestIndex_offsetPost] = min(differences4);

                clear tmp
                tmp=pupil_area_zscore(nearestIndex_onset:nearestIndex_offset);
                pupil_trials_stim(istim,isf, itrial,1:length(tmp))=tmp;

                clear tmp2
                tmp2=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_onset);
                pupil_trials_blanks(istim,isf, itrial,1:length(tmp2))=tmp2;

                clear tmp3
                tmp3=pupil_area_zscore(nearestIndex_offset:nearestIndex_offsetPost);
                pupil_trials_poststim(istim,isf, itrial,1:length(tmp3))=tmp3;

                clear tmp4
                tmp4=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_offset);
                pupil_trials_epochs(istim,isf, itrial,1:length(tmp4))=tmp4;

                clear tmp5
                tmp5=pupil_area_zscore(nearestIndex_onsetPre:nearestIndex_offsetPost);
                pupil_trials_entire_epochs(istim,isf, itrial,1:length(tmp5))=tmp5;

            end  
       end
    
end
pupil_epochs.pupil_trials_entire_epochs=pupil_trials_entire_epochs
pupil_epochs.pupil_trials_blanks=pupil_trials_blanks
pupil_epochs.pupil_trials_stim=pupil_trials_stim

pupil_epochs.pupil_trials_poststim=pupil_trials_poststim
pupil_epochs.pupil_trials_epochs=pupil_trials_epochs

pupil_epochs.stimulus=stim_directions
pupil_epochs.dirs_stimulus=stim_directions
pupil_epochs.sfs_stimulus=stim_sfs
pupil_epochs.ntrials=ntrials
pupil_epochs.prestim_duration=prestim_duration
pupil_epochs.poststim_duration=poststim_duration
pupil_epochs.stim_duration=stim_duration
pupil_epochs.type='fixed sequence'
pupil_epochs.stimulus_duration_in_sec='5s'

pupil_mat_name = [filepathanalysis, '\pupil_epochs.mat'];
save(pupil_mat_name, 'pupil_epochs');


