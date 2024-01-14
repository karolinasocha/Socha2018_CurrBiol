function expt=doLoadStimLogs3(expt,varargin)
%DOLOADSTIMLOGS
% EXPT=DOLOADSTIMLOGS(EXPT)

% 140923 SK, elseif to stimlog added

defaultopts = {'Overwrite',false,'FixStimLogs',false};

options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1};
end

fnStimLog = fullfile(expt.dirs.analrootpn,'exptstimlog.mat');

if exist(fnStimLog) & ~options.Overwrite
    fprintf('Loading cached logs file %s\n',fnStimLog);
    expt = load(fnStimLog);
    
    return;
end

% reads and parse stimulatio logs
if isfield(expt.filenames,'presentationlog')
    fprintf(1,'Reading presentation log\n');
    stimlog.raw = readPresentationLog(expt.filenames.presentationlog);
    
    fprintf(1,'Parsing presentation log\n');
    expt.stimlog.parsed = parsePresentationLog3(stimlog.raw); %to fix skipping frames
%     expt.stimlog.parsed = parsePresentationLog2(stimlog.raw);
%     expt.stimlog.fixed = fixInterrupedLog(expt.stimlog.parsed);
    expt.stimlog.fixed = expt.stimlog.parsed;
    
    expt.info = expt.stimlog.fixed.info;
    expt.nFramesPerTrial = expt.info.nStim*(expt.info.nOffPulses+expt.info.nOnPulses);
    expt.nFramesTotal = expt.info.nTrials*expt.nFramesPerTrial;

    [expt.frames.stims,expt.frames.blanks, expt.frames.epochs]=getPresentationEpochs(expt.stimlog.fixed);
    
    expt.dirs.stimframes = fullfile(expt.dirs.stimulation,'bitmaps',expt.info.animal,expt.info.expt);
    
    expt.frameRate = ((expt.stimlog.parsed.FrameTimes(end)-expt.stimlog.parsed.FrameTimes(1))/(length(expt.stimlog.parsed.FrameTimes)-1)/1e4).^-1;
    
% % elseif isfield(expt.filenames,'stimlog') % SK 
% %     fprintf(1,'Reading presentation log\n');
% %     stimlog.raw = readPresentationLog(fullfile(expt.dirs.analysis,expt.animal,expt.name,expt.filenames.stimlog));
% %     
% %     fprintf(1,'Parsing presentation log\n');
% %     expt.stimlog.parsed = parsePresentationLog2(stimlog.raw);
% % %     expt.stimlog.fixed = fixInterrupedLog(expt.stimlog.parsed);
% %     expt.stimlog.fixed = expt.stimlog.parsed;
% %     
% %     expt.info = expt.stimlog.fixed.info;
% %     expt.nFramesPerTrial = expt.info.nStim*(expt.info.nOffPulses+expt.info.nOnPulses);
% %     expt.nFramesTotal = expt.info.nTrials*expt.nFramesPerTrial;
% % 
% %     [expt.frames.stims,expt.frames.blanks, expt.frames.epochs]=getPresentationEpochs(expt.stimlog.fixed);
% %     
% %     expt.dirs.stimframes = fullfile(expt.dirs.stimulation,'bitmaps',expt.info.animal,expt.info.expt);
% %     
% %     expt.frameRate = ((expt.stimlog.parsed.FrameTimes(end)-expt.stimlog.parsed.FrameTimes(1))/(length(expt.stimlog.parsed.FrameTimes)-1)/1e4).^-1;
    
else
    fprintf('Stimulation log not found!\n');
    expt.nFramesTotal = inf;
end

% reads stimulation protocol
if isfield(expt,'info')
    
    expt.dirs.prots = fullfile(expt.dirs.stimulation,'prots',expt.info.animal);
    temp = fullfile(expt.dirs.prots,[expt.info.expt '.prot']);
    
    if exist(temp)==2
        expt.filenames.prot = temp;
        expt.prot = readprot(expt.filenames.prot);
    end
    
else
    fprintf('Protocol not found!\n');    
end

if isfield(expt.dirs,'stimulation') && isfield(expt,'info')

    temp = fullfile(expt.dirs.stimulation,'bitmaps',expt.animal,expt.info.expt);
    
    if exist(temp)==7
        expt.dirs.stims= temp;
    end
else
    fprintf('Stimulus frames not found!\n');    
end

if ~isfield(expt,'desc')
    expt.desc = expt.id;
end

if isfield(expt,'prot') & isstruct(expt.prot)
    expt.desc = [expt.desc ' ' expt.prot.StimulusType];
end

if isfield(expt,'info') & isstruct(expt.info)
    expt.desc = [expt.desc ' ' expt.info.animal '/' expt.info.expt];
end

if options.FixStimLogs
    expt = doFixStimLogs(expt);
end

fprintf('Saving log cache %s \n',fnStimLog);
save(fnStimLog,'-struct','expt')

return;