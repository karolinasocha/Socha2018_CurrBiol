function [av_eye_epochs,trials_eye_epochs, av_eye_stims, trials_eye_stims]=calculate_average_pupil(eye_diam_entire_trace, epochs, stim_epochs)

try
[av_eye_epochs trials_eye_epochs]=tcEpochAverage2(eye_diam_entire_trace,epochs);% time missing 
end
try
[av_eye_epochs trials_eye_epochs]=tcEpochAverage2(eye_diam_entire_trace',epochs);% time missing 
end
try
[av_eye_stims trials_eye_stims]=tcEpochAverage2(eye_diam_entire_trace,stim_epochs);% time missing 
end
try
[av_eye_stims trials_eye_stims]=tcEpochAverage2(eye_diam_entire_trace',stim_epochs);% time missing 
end
