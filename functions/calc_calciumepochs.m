function [ddat ddat_stim]=calc_calciumepochs(expt2,calcium_trace)

ntrials=expt2.info.nTrials;
nstims=expt2.info.nStim;
frameRate=expt2.frameRate;
cframespre=round(2*frameRate); % 2 sec pre
stimduration=expt2.info.nOnPulses;
cframespost = stimduration+round(3*frameRate); % 3 sec
length_epoch=cframespost+cframespre;
n_boutons=size(calcium_trace,2);
ddat = nan(n_boutons,length_epoch,ntrials,nstims);
ddat_stim=nan(n_boutons,stimduration,ntrials,nstims);
size(ddat);

X=calcium_trace;
%Z_velocity=lowpass(X);
for iStim = 1:nstims-1
    for iTrial = 1:ntrials              
          camidx=expt2.frames.stims{iTrial,iStim}(1);
          camidx_end=expt2.frames.stims{iTrial,iStim}(end);
            for ibouton=1:n_boutons
                ddat(ibouton,:,iTrial,iStim) = X(1-cframespre + camidx:camidx+cframespost,ibouton);  
                ddat_stim(ibouton,1:length(camidx:camidx_end),iTrial,iStim)=X(camidx:camidx_end,ibouton);
            end
    end
end

return