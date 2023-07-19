
function [trials_eye_epochs_reordered,trials_eye_stims_reordered]=calculate_same_stimulus_order(expt2,stimulus,trials_eye_epochs, trials_eye_stims)


    clear correct_value
    clear order
    frameRate=round(expt2.frameRate);
    stim_epochs_length=length(trials_eye_stims);
    all_epochs_length=length(trials_eye_epochs);

    stims_value=stimulus(1:end-1);
    correct_value=[0,30,60,90,120,150,180,210,240,270,300,330];

    for i=1:12;
    order(i)= find(stims_value==correct_value(i));
    end

    if length(size(trials_eye_stims))==4
        trials_eye_stims_reordered=trials_eye_stims(:,:,:,order);
        trials_eye_epochs_reordered=trials_eye_epochs(:,:,:,order);

    elseif length(size(trials_eye_stims))==2
        trials_eye_stims_reordered=trials_eye_stims(:,order);
        trials_eye_epochs_reordered=trials_eye_epochs(:,order);
        
    else
        disp('Wrong size of matrix')
end