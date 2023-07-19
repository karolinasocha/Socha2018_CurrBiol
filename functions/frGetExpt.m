function expt = frGetExpt(newdir,mode)
%FRGETEXPT Creates experiment structure
% EXPT = FRGETEXPT(NEWDIR) where NEWDIR is string with 
% format $ANIMAL\$EXPTNAME

% 10/01 id uses windows directory format

% 140911 SK, ica mask path added
% 140904 SK, rec path changed, delete option with reg
% 131004 SK, eye/behav movie paths added
% 131107 SK, aligned mask path added
% 131202 SK, presentation path corrected, undone
% 140210 SK, windowcam paths added
% 140222 SK, astro paths added

[animal,exptname]=fileparts(newdir);

if nargin < 2
    mode = '2photon';
end

expt.animal = animal;
expt.name = exptname;

expt.id = sprintf('%s\\%s',expt.animal,expt.name);

expt.dirs = frGetDirs;

%% image root
expt.dirs.rawroot1p = fullfile(expt.dirs.widefield,expt.animal);
expt.dirs.rawrootpn = fullfile(expt.dirs.raw,expt.animal,expt.name);
expt.dirs.recrootpn = fullfile(expt.dirs.reconstructed,expt.animal,expt.name);
expt.dirs.regrootpn = fullfile(expt.dirs.registered,expt.animal,expt.name);

if exist(expt.dirs.rawrootpn)~=7 
    expt.dirs.rawrootpn = fullfile(expt.dirs.backup,'2photon\raw',expt.animal,expt.name);
end;

if exist(expt.dirs.rawrootpn)~=7 
    expt.dirs.rawrootpn = fullfile('Q:\data','2photon\raw',expt.animal,expt.name);
end;

if exist(expt.dirs.recrootpn)~=7 
    disp('don''t touch ! if this breaks your code ask vincent');
%     expt.dirs.recrootpn = fullfile(expt.dirs.backup,'2photon\rec',expt.animal,expt.name); % nerffs version
    expt.dirs.recrootpn = fullfile(expt.dirs.reconstructed,expt.animal,expt.name); % for deleting rec files on bkrunch
end;

if exist(expt.dirs.regrootpn)~=7 
    disp('don''t touch ! if this breaks your code ask vincent');
%     expt.dirs.regrootpn = fullfile(expt.dirs.backup,'2photon\reg',expt.animal,expt.name);
%     disp('Steffen changed this for temp use of brkunch for registering data')
    expt.dirs.regrootpn = fullfile(expt.dirs.registered,expt.animal,expt.name); 
    if ~exist(expt.dirs.regrootpn) %KM 150508
        expt.dirs.regrootpn = fullfile(expt.dirs.backup,'2photon\reg',expt.animal,expt.name); %KM 150508
        disp('If you see this, your file is not in H: It loaded from nerffs01 --KM')
    end
end;

% %% raw data directories

% scanlog = frParseScanLog(expt.dirs.rawrootpn); 
% 
% if isempty(scanlog)   
%     pn = fullfile(expt.dirs.data,expt.animal,expt.name);
%     if exist(pn)==7
%         expt.dirs.dadsrootpn = pn;
%         % read scan log if exists
%         scanlog = frParseScanLog(expt.dirs.dadsrootpn);    
%     end
% end
% 
% expt.scanlog = scanlog;

% if ~isempty(expt.scanlog)
%     expt.desc = [expt.scanlog.start_date ' ' expt.id];
% else
%     expt.desc = expt.id;
% end

%% raw directories

%% green / red image directories
%expt.dirs.reggreenpn = expt.dirs.regrootpn; % temporary added KS
if ~isdir(expt.dirs.rawrootpn)
    str = sprintf('%s not found',expt.dirs.rawrootpn);
    warning(str);
end;

if ~isdir(expt.dirs.recrootpn)
    str = sprintf('%s not found',expt.dirs.recrootpn);
    warning(str);
end;

if isdir(expt.dirs.rawroot1p)
    str = sprintf('%s not found',expt.dirs.rawroot1p);
    warning(str);
end;

list = dir(expt.dirs.rawrootpn);
list2 = dir(expt.dirs.recrootpn);
list3 = dir(expt.dirs.regrootpn);
list4 = dir(expt.dirs.rawroot1p);
    
if sum(strmatch('.tif',{list.name}))
    expt.dirs.rawgreenpn = expt.dirs.rawrootpn;
    expt.dirs.rawredpn = '';
    expt.dirs.recgreenpn = expt.dirs.recrootpn;
    expt.dirs.reggreenpn = expt.dirs.regrootpn;
elseif sum(strmatch('.tif',{list2.name}))
    expt.dirs.rawredpn = '';
    expt.dirs.recgreenpn = expt.dirs.recrootpn;
    expt.dirs.reggreenpn = expt.dirs.regrootpn;
elseif sum(strmatch('.tif',{list3.name}))
    expt.dirs.reggreenpn = expt.dirs.regrootpn;    
elseif sum(strmatch('.tif',{list4.name}))
    expt.dirs.rawredpn = '';
    expt.dirs.recgreenpn = '';
    expt.dirs.reggreenpn = '';
else    
    warning('no tif files found');
end;

%% analysis directories
if nargin < 3
    expt.dirs.analrootpn = fullfile(expt.dirs.analysis,expt.animal,expt.name);
else
    expt.dirs.analrootpn = fullfile(expt.dirs.analysis,expt.animal,expt.name,expt.planeid);
end;

if exist(expt.dirs.analrootpn)~=7
    [success]=mkdir(expt.dirs.analrootpn);
    if ~success
        warning(sprintf('Could not create directory: %s',expt.dirs.analrootpn))
    end
end;

%% analysis directories
if nargin < 3
    expt.dirs.figrootpn = fullfile(expt.dirs.figures,expt.animal,expt.name);
else  % case with multiple planes
    expt.dirs.figrootpn = fullfile(expt.dirs.figures,expt.animal,expt.name,expt.planeid);
end;

if exist(expt.dirs.figrootpn)~=7
    [success]=mkdir(expt.dirs.figrootpn);
    if ~success
        warning(sprintf('Could not create directory: %s',expt.dirs.analrootpn))
    end;
end;

%% standard filenames for derived data

% uncomment below if actually needed

expt.filenames.registered = [exptname,'_green_reg.tif'];
% expt.filenames.pca = [exptname,'_green_pca.tif'];
% expt.filenames.smoothed = [exptname,'_green_reg_sm.tif'];
% expt.filenames.registeredred = [exptname,'_red_reg.tif'];

expt.filenames.stack_downsampled = 'stack_downsampled.tif';
expt.filenames.shifts = [exptname,'_shifts.mat'];
expt.filenames.masks = 'masks_neurons.mat';
expt.filenames.masks_KS = 'masks.mat'; %KS
expt.filenames.masksZoom = 'masks_neuronsZoom.mat';
expt.filenames.masks_auto = 'masks_neurons_auto.mat';
expt.filenames.masks_aligned =  'masks_neurons_aligned.mat'; % SK
expt.filenames.masks_active =  'masks_neurons_active.mat'; % SK
expt.filenames.masks_astros =  'masks_astros.mat'; % SK
expt.filenames.masks_domains =  'masks_domains.mat'; % SK
expt.filenames.masks_ica =  'masks_ica.mat'; % SK
expt.filenames.timecourses = 'timecourses.mat';
expt.filenames.timecoursespca = 'timecourses_pca.mat';
expt.filenames.timecourses_astros = 'timecourses_astros.mat'; % SK
expt.filenames.timecourses_domains = 'timecourses_domains.mat'; % SK
expt.filenames.pca_usv = 'pca_usv.mat';
expt.filenames.imgs = 'imgs.mat';
expt.filenames.deconv= 'deconv.mat';
expt.filenames.revcorr = 'revcorr.mat';
expt.filenames.stimlog = 'stimlog.mat';
expt.filenames.maps='maps.mat'; %KS
expt.filenames.tcs_KS='tcs.mat'; %KS
expt.filenames.roiwindow='roiwindow.mat'; % added 150213

%% visual stimulation directories
temp = fullfile(expt.dirs.presentation,expt.animal);

% presentation log directory
if exist(temp)==7 
    expt.dirs.presentation = temp;
end

% presentation log file
% temp = fullfile(expt.dirs.presentation,expt.animal,[expt.name,'.log']); % SK expt.animal added
temp = fullfile(expt.dirs.presentation,[expt.name,'.log']);
if exist(temp)==2
    expt.filenames.presentationlog = temp;
else
    disp('warning: could not find presentation log');
end

%% cam data, SK added paths to cam data
% eyecam
tmp = fullfile(expt.dirs.eyes,expt.animal);
if exist(tmp) == 7
    expt.dirs.eyes = tmp;
end
tmp = fullfile(expt.dirs.eyes, [expt.name,'.seq']);
if exist(tmp) == 2
    expt.filenames.eyemovie = tmp;
else
%     expt.filenames.eyemovie = '';
    tmp2 = fullfile(expt.dirs.backup,'eyecam',expt.animal,[expt.name '.seq']);
    if exist(tmp2) == 2
        expt.filenames.eyemovie = tmp2;
        disp('warning: setting backup eyecam path');
    else
        expt.filenames.eyemovie = [];
        disp('warning: could not find eye movie');
    end
end

% facecam
tmp = fullfile(expt.dirs.behavior,expt.animal);
if exist(tmp) == 7
    expt.dirs.behavior = tmp;
end
tmp = fullfile(expt.dirs.behavior, [expt.name,'.seq']);
if exist(tmp) == 2
    expt.filenames.behaviormovie = tmp;
else
%     expt.filenames.behaviormovie = '';
    tmp2 = fullfile(expt.dirs.backup,'facecam',expt.animal,[expt.name '.seq']);
    if exist(tmp2) == 2
        expt.filenames.behaviormovie = tmp2;
        disp('warning: setting backup facecam path');
    else
        expt.filenames.behaviormovie = [];
        disp('warning: could not find behavior movie'); 
    end
end

% windowcam
tmp = fullfile(expt.dirs.window,expt.animal);
if exist(tmp) == 7
    expt.dirs.behavior = tmp;
end
tmp = fullfile(expt.dirs.window, [expt.name,'.seq']);
if exist(tmp) == 2
    expt.filenames.windowmovie = tmp;
else
%     expt.filenames.windowmovie = '';
    tmp2 = fullfile(expt.dirs.backup,'windowcam',expt.animal,[expt.name '.seq']);
    if exist(tmp2) == 2
        expt.filenames.windowmovie = tmp2;
        disp('warning: setting backup windowcam path');
    else
        expt.filenames.windowmovie = [];
        disp('warning: could not find window movie'); 
    end
end

return;
