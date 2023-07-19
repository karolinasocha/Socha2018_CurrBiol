function dirs = frGetDirs
% FRGETDIRS Returns directory structure based on user name and host name.
%
% by Vincent Bonin based on directories.m by Mark Histed
% $id$
%
% 140122 SK, H/Q drive dirs prep (commented)
% 140904 SK, cleanup
% 131002 SK, dirs added
% 140210 SK, dirs added
% 140430 SK, dropbox path edited

%% KS commented this section since it crashes function 04-07-2023
% global USER_PROFILE;
% 
% if isempty(USER_PROFILE)
%     user = username;
% else
%     user = USER_PROFILE;
% end;

%%
hostname=getenv('COMPUTERNAME'); %210709, replace outdated hostname.m
user=getenv('username');

switch lower(hostname)
    case 'bkrunch'
        switch user
            case {'Administrator',...
                    'cagatay',...
                    'daniel',...
                    'dun',...
                    'alexander',...
                    'jessica',...
                    'kaan',...
                    'karolina',...
                    'katrien',...
                    'michal',...
                    'mouselab',...
                    'vincent',...
                    'wclee',...
                    'xuhan',...
                    'filip',...
                    }
                % dirs.root = 'h:\users\vincent';
                dirs.svn = 'g:\mousebox\code\mouselab';
                dirs.analysis = 'g:\mousebox\analysis';
                dirs.figures = 'g:\mousebox\figures'; % SK
                dirs.scripts = ''; % SK
                dirs.data = 'h:\data';
                dirs.local = 'h:\local';
                dirs.data2 = 'q:\data';
                dirs.local2 = 'q:\local';
                % dirs.stimulation = fullfile(dirs.data,'vstim'); % XH
                dirs.stimulation = fullfile('h:\data','vstim'); % XH, in case of conflicting with q:\
                dirs.presentation = fullfile(dirs.data,'presentation');
                dirs.raw = fullfile(dirs.data,'2photon\raw');
                dirs.reconstructed = fullfile(dirs.data,'2photon\rec');
                dirs.registered = fullfile(dirs.data,'2photon\reg');
                dirs.widefield = fullfile(dirs.data,'1photon\raw'); % SK
                dirs.eyes = fullfile(dirs.data,'eyecam'); % SK
                dirs.behavior = fullfile(dirs.data,'facecam'); % SK
                dirs.window = fullfile(dirs.data,'windowcam'); % SK
                % dirs.pca = 'h:\data\pca';
                % dirs.random = 'h:\data\vis_stim_protocols';
                % dirs.texture = fullfile(dirs.svn,'thirdparty\textureSynth'); % KS
                % dirs.pyr = fullfile(dirs.svn,'thirdparty\thirdparty\matlabPyrTools');% KS
                dirs.backup = '\\nerffs01\mouselab\data';
                
            case 'steffen'
                dirs.root = 'h:\local\users\steffen';
                dirs.svn = 'g:\mousebox\code\mouselab';
                dirs.analysis = 'g:\mousebox\analysis';
                dirs.figures = 'g:\mousebox\figures';
                dirs.scripts = '';
                dirs.backup = '\\nerffs01\mouselab\data';
                dirs.data = 'h:\data';
                dirs.local = 'h:\local';
                dirs.stimulation = fullfile(dirs.backup,'vstim');
                dirs.presentation = fullfile(dirs.backup,'presentation');
                dirs.raw = fullfile(dirs.backup,'2photon\raw');
                dirs.reconstructed = fullfile(dirs.data,'2photon\rec');
                dirs.registered = fullfile(dirs.backup,'2photon\reg');
                dirs.widefield = fullfile(dirs.backup,'1photon\raw');
                dirs.eyes = fullfile(dirs.backup,'eyecam');
                dirs.behavior = fullfile(dirs.backup,'facecam');
                dirs.window = fullfile(dirs.backup,'windowcam');
                
                
            case {'ben'} % BV20160105: created separate case so I can play with folder and stuff, without breaking things for other users
                dirs.root = 'h:\local\users\ben';
                dirs.svn = 'g:\mousebox\code\mouselab';
                dirs.analysis = 'g:\mousebox\analysis';
                dirs.figures = 'g:\mousebox\figures';
                dirs.scripts = '';
                dirs.backup = '\\nerffs01\mouselab\data';
                dirs.data = 'q:\data';
                dirs.local = 'q:\local';
                dirs.stimulation = fullfile(dirs.data,'vstim');
                dirs.presentation = fullfile(dirs.data,'presentation');
                dirs.raw = fullfile(dirs.data,'2photon\raw');
                dirs.reconstructed = fullfile(dirs.data,'2photon\rec');
                dirs.registered = fullfile(dirs.data,'2photon\reg');
                dirs.widefield = fullfile(dirs.data,'1photon\raw');
                dirs.eyes = fullfile(dirs.data,'eyecam');
                dirs.behavior = fullfile(dirs.data,'facecam');
                dirs.window = fullfile(dirs.data,'windowcam');
                
            otherwise
                error('unknown user name');
        end
        
    otherwise
        error('unknown host name');
end

return;
