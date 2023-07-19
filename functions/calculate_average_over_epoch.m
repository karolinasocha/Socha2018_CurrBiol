function [av_eye_stims_reordered_onset, trials_eye_stims_reordered_onset,...
    av_eye_stims_reordered_offset, trials_eye_stims_reordered_offset,...
    av_eye_stims_reordered_middle, trials_eye_stims_reordered_middle,...
    av_eye_stims_reordered_entire, trials_eye_stims_reordered_entire]= calculate_average_over_epoch(expt2, av_eye_stims_reordered, trials_eye_stims_reordered)

frameRate=expt2.frameRate;
frame_starts=round(0.5*frameRate)
frame_ends=round(1*frameRate)

% onsets average between 0.5sec and 1sec. when stimulus starts
av_eye_stims_reordered_onset=squeeze(nanmean(av_eye_stims_reordered(frame_starts:frame_ends, :),1));
trials_eye_stims_reordered_onset=squeeze(nanmean(trials_eye_stims_reordered(frame_starts:frame_ends,:,:,:),1));


% offset averaged at the last 0.5 sec. of stimulus duration
av_eye_stims_reordered_offset=squeeze(nanmean(av_eye_stims_reordered(end-frame_starts:end, :),1));
trials_eye_stims_reordered_offset=squeeze(nanmean(trials_eye_stims_reordered(end-frame_starts:end,:,:,:),1));

% middle value between 1.5sec - 4.5sec. of stimulus duration
av_eye_stims_reordered_middle=squeeze(nanmean(av_eye_stims_reordered(frame_ends+1:end-frame_starts-1, :),1));
trials_eye_stims_reordered_middle=squeeze(nanmean(trials_eye_stims_reordered(frame_ends+1:end-frame_starts-1,:,:,:),1));

% entire epoch
av_eye_stims_reordered_entire=squeeze(nanmean(av_eye_stims_reordered(:,:),1));
trials_eye_stims_reordered_entire=squeeze(nanmean(trials_eye_stims_reordered(:,:,:,:),1));

% %not sure of useful
% data_pupil.av_eye_stims_reordered_onset=av_eye_stims_reordered_onset;
% data_pupil.trials_eye_stims_reordered_onset=trials_eye_stims_reordered_onset;
% 
% data_pupil.av_eye_stims_reordered_offset=av_eye_stims_reordered_offset;
% data_pupil.trials_eye_stims_reordered_offset=trials_eye_stims_reordered_offset;
% data_pupil.av_eye_stims_reordered_middle=av_eye_stims_reordered_middle;
% data_pupil.trials_eye_stims_reordered_middle=trials_eye_stims_reordered_middle;
% 
% data_pupil.trials_eye_stims_reordered_entire=trials_eye_stims_reordered_entire;
% data_pupil.av_eye_stims_reordered_entire=av_eye_stims_reordered_entire;

end
