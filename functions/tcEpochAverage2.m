function [av,trials, sd] = tcEpochAverage2(timeCourses,epochs, ratio)
%TCTRIALAVERAGE
% [AV,TRIALS, SD]= TCEPOCHAVERAGE2(TIMECOURSES,EPOCHS)
% [AV,TRIALS, SD]= TCEPOCHAVERAGE2(TIMECOURSES,EPOCHS, RATIO)

if nargin < 3
    ratio = 1;
end;

eLen = getEpochsLen(epochs);

[nSamples,nCells] = size(timeCourses);

[nTrials,nStim]=size(epochs);

trials = zeros(eLen/ratio,nCells,nTrials,nStim);

for iTrial = 1:nTrials
    for iStim = 1:nStim
        
        ind = epochs{iTrial,iStim}(1:eLen);   
     
        temp =  timeCourses(ind,:);
        
        if ratio > 1
            temp = tcDecimate(temp,ratio);
        end;

       	trials(:,:,iTrial,iStim) = temp;
        
    end
end

av = squeeze(mean(trials,3));
sd = squeeze(std(trials,[],3));

return;