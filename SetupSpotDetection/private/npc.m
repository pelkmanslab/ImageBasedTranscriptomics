function localPath = npc(networkPath)
%NPC is a network path conversion function mapping remote path to its local
%counterpart.
%
% This version contian a patch by Thomas:
% x Does not give constant misleading warnings
% x Uses same optional settings file as npc
%
% Network Path Conversion function is used to wrap around inline string
% constants with full path to a file. Input is always a path as it is valid
% for the remote side (e.g. cluster, etc.) and it should be substituted
% with a valid local path (usable on the machine which calls this function)
% based on settings file called 'npc.local.mat'
%
% IMPORTANT! 'npc.local.mat' should be never committed into the code
% repository.
%
% Settings
%
%     settings = struct();
%     settings(1).pattern = '/BIOL/imsb/fs2/bio3/bio3';
%     settings(1).replace = '/share/nas/ethz-share2';
%     settings(2).pattern = '/BIOL/imsb/fs3/bio3/bio3';
%     settings(2).replace = '/share/nas/ethz-share3';
%     settings(3).pattern = '/BIOL/sonas/biol_uzh_pelkmans_s4';
%     settings(3).replace = '/share/nas/ethz-share4';
%     settings(4).pattern = '/BIOL/sonas/biol_uzh_pelkmans_s6';
%     settings(4).replace = '/share/nas/ethz-share6';
%     settings(5).pattern = '/BIOL/sonas/biol_uzh_pelkmans_s7';
%     settings(5).replace = '/share/nas/ethz-share7';
%     filename = [os.path.dirname(which('npc')) filesep 'npc.local'];
%
% Then to save as .MAT
%
%     save([filename '.mat'],'settings');
%
% or as
%
%     string.write([filename '.json'], savejson('', settings));
%
%
% @author: Berend Snijder <berend.snijder@imls.uzh.ch>
% @author: Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
%
% Lucas Pelkmans Lab http://www.pelkmanslab.org/
%
% See also FULLFILE, FILESEP, OS.PATH.EXISTS


if nargin==0
    % the formats we want to be able to convert
    
    %     networkPath = '\\nas-biol-imsb-1\share-3-$\Data\Users\Berend';
    %     networkPath = 'X:\Data\Users\Berend';
    %     networkPath = '/Volumes/share-3-$/Data/Users/Pauli/';
    %     networkPath = '\BIOL\imsb\fs3\bio3\bio3\Data\Users\Prisca\090203_Mz_Tf_EEA1\090203_Mz_Tf_EEA1_CP395-1ad\TIFF\TIFF_B05_s20_w1608AC8BF-29E4-498C-AA60-BBCFE519CB35.png';
    %     networkPath = 'http://www.ibrain.ethz.ch/share-3/Data/Users/Prisca/Safia/iBRAIN/38exp/38exp-R3/BATCH';
    %     networkPath = '\\nas-biol-imsb-1.d.ethz.ch\share-3-$\Data\Users\Prisca\090203_Mz_Tf_EEA1_vesicles\090203_Mz_Tf_EEA1_CP392-1ad\SEGMENTATION\';
    %     networkPath = '/IMSB/images/thaminys/iBrain/45exp/45expR1/BATCH/';
    %     networkPath = '/IMSB/images/thaminys/iBrain/21exp/magnification_10x/time_0hour/BATCH/Batch_data.mat';
    networkPath = '/BIOL/imsb/fs2/bio3/bio3/Data/Users/Berend';
    %     networkPath = os.path.dirname(mfilename('fullpath'));
end

% Map path using local settings (smart way of doing this is to use
% persitent variables - see invoked functions for details).
%
% save('.../path/npc.local.mat', 'settings') where settings is a struct
% array with 'pattern' and 'replace' keys.
%
localSettings = getLocalSettings();
% Has local settings? Use the map the path. Note: remote sites like Brutus
% never have any local settings.
if ~isempty(localSettings)
    localPath = mapNetworkPath(networkPath, localSettings);
    try
        if strcmp(localPath,networkPath)
            % Also locally use heuristics, if no change
            localPath = fallback(networkPath);
        end
    catch exception
        warning('Error during fallback: ignore fallback due to the error: %s', exception.message);
    end
else
    % Fallback to original heuristics.
    localPath = fallback(networkPath);
end

if ~os.path.exists(localPath)
    warning('npc: Returning a non-existing local path "%s"', localPath);
end

end



function settings = getLocalSettings()
persistent localSettings;
if isempty(localSettings)
    localSettings = [];
    [pathstr, name, ext] = fileparts(mfilename('fullpath'));
    matSettingsFile = [pathstr filesep 'npc.local.mat'];            % use same setting file as npc
    jsonSettingsFile = [pathstr filesep 'npc.local.json'];
    if os.path.exists(matSettingsFile)
        data = load(matSettingsFile);
        localSettings = data.settings;
    elseif os.path.exists(jsonSettingsFile)
        if ~exist('loadjson')
            % Please install JSONlab package from
            % http://www.mathworks.com/matlabcentral/fileexchange/33381
            settings = [];
            return
        end
        localSettings = loadjson(jsonSettingsFile);
    end
end
settings = localSettings;
end


function convPath = mapNetworkPath(thePath, settings)
convPath = thePath;
for index = 1:numel(settings)
    if strfind(thePath, settings(index).pattern)
        convPath =  strrep(...
            thePath, settings(index).pattern, settings(index).replace);
        return
    end
end
end


function strConvPath = fallback(strRootPath)
% by default, the output is the same as the input, only if the input could
% be succesfully converted do we set the conversion result as output
strConvPath = strRootPath;

% detect format of input path, and store value in intPathFormat, possible
% values are:
% intPathFormat = 1 --> PC (mapped),
% intPathFormat = 2 --> PC,
% intPathFormat = 3 --> Mac
% intPathFormat = 4 --> UNIX/Brutus
% intPathFormat = 5 --> HTTP
% intPathFormat = 6 --> UNIX/Brutus for Safia /IMSB/images

intPathFormat = 0;

% in case we're dealing with an old Hreidar path, replace the initial
% mapping to the NAS... (can only have been share-2-$, share-3-$ didnt
% exist in Hreidar days)
if strncmpi(strRootPath,'/hreidar/extern/bio3/',21)
    strRootPath = strcat('/BIOL/imsb/fs2/bio3/bio3/',strRootPath(22:end));
end


% see what type of input we're dealing with
% mapped PC path references (case insensitive)
if ~isempty(regexpi(strRootPath,'^[E-Z]:\\'))
    intPathFormat = 1; % PC format mapped
    strShareRegexp = '^[E-Z]:\\'; % --> this helps to catch the PathBase, but not the share number...
    
    % if it's a mapped windows path, see if the current system mapping is
    % valid for this drive, if so, this is the most likely case and return
    strCurrentlyMappedPath = replaceDRIVEwithUNC(strRootPath);
    if fileattrib(strCurrentlyMappedPath)
        strConvPath = strCurrentlyMappedPath;
        return
    end
    
    %     disp('input path is mapped PC')
    
    % typical PC path
elseif ~isempty(regexpi(strRootPath,'^\\\\nas-biol-imsb-1'))
    intPathFormat = 2; % PC format
    strShareRegexp = '^\\\\nas-biol-imsb-1.*\\share-(\d)-\$\\';
    %     disp('input path is PC')
    
elseif ~isempty(regexpi(strRootPath,'^\\\\nas-unizh-imsb1'))
    intPathFormat = 2; % PC format
    strShareRegexp = '^\\\\nas-unizh-imsb1.*\\share-(\d)-\$\\';
    %     disp('input path is PC')
    
elseif ~isempty(regexpi(strRootPath,'^\\\\nas21nwg01.ethz.ch'))
    intPathFormat = 2; % PC format
    strShareRegexp = '^\\\\nas21nwg01.ethz.ch\\biol_uzh_pelkmans_s(\d)\\';
    %     disp('input path is PC')
    
    
    % typical Mac path references (not SONAS)
elseif ~isempty(regexp(strRootPath,'^/Volumes/share-'))
    intPathFormat = 3; % Mac format
    strShareRegexp = '^/Volumes/share-(\d)-\$/';
    %     disp('input path is Mac')
    
    % typical Mac path references (SONAS)
elseif ~isempty(regexp(strRootPath,'^/Volumes/biol_uzh_pelkmans_s'))
    intPathFormat = 3; % Mac format
    strShareRegexp = '^/Volumes/biol_uzh_pelkmans_s(\d)/';
    %     disp('input path is Mac')
    
    % alternative Mac path references
    %(please see WIKI for shell script to mount shares to defined mount points in the home directory on mac)
elseif ~isempty(regexp(strRootPath,'.*/shares/ethz-share'))
    intPathFormat = 3; % Mac format
    strShareRegexp = '.*shares/ethz-share(\d)/'; % create user specific path
    %     disp('input path is Mac')
    
    % typical Brutus cluster path references (not SONAS)
elseif ~isempty(regexp(strRootPath,'^/BIOL/imsb'))
    intPathFormat = 4; % Brutus format
    strShareRegexp = '^/BIOL/imsb/fs(\d)/bio3/bio3/';
    %     disp('input path is UNIX/Brutus')
    
    % weird brutus paths... (not SONAS)
elseif ~isempty(regexp(strRootPath,'^\\BIOL\\imsb'))
    strRootPath = strrep(strRootPath,'\','/');
    intPathFormat = 4; % Brutus format
    strShareRegexp = '^/BIOL/imsb/fs(\d)/bio3/bio3/';
    %     disp('input path is UNIX/Brutus')
    
    % typical Brutus cluster path references (SONAS)
elseif ~isempty(regexp(strRootPath,'^/BIOL/sonas/biol_uzh_pelkmans_s'))
    intPathFormat = 4; % Brutus format
    strShareRegexp = '^/BIOL/sonas/biol_uzh_pelkmans_s(\d)/';
    %     disp('input path is UNIX/Brutus')
    
    % weird brutus paths (SONAS)...
elseif ~isempty(regexp(strRootPath,'^\\BIOL\\sonas\\biol_uzh_pelkmans_s'))
    strRootPath = strrep(strRootPath,'\','/');
    intPathFormat = 4; % Brutus format
    strShareRegexp = '^/BIOL/sonas/biol_uzh_pelkmans_s(\d)/';
    %     disp('input path is UNIX/Brutus')
    
    % webpage paths: original website
elseif ~isempty(regexp(strRootPath,'^http://www\.ibrain\.ethz\.ch/share-'))
    intPathFormat = 5; % HTTP website format
    strShareRegexp = '^http://www\.ibrain\.ethz\.ch/share-(\d)/';
    
    % webpage paths: pelkmanslab.org
elseif ~isempty(regexp(strRootPath,'^https://ibrain\.pelkmanslab\.org/share-'))
    intPathFormat = 5; % HTTP website format
    strShareRegexp = '^https://ibrain\.pelkmanslab\.org/share-(\d)/';
    
    
elseif ~isempty(regexp(strRootPath,'^/IMSB/images', 'once'))
    % SAFIAs disk is a special case. hanlde it here.
    if ispc
        strNewPath = strrep(strRootPath,'/IMSB/images/','\\nas-biol-ibt-1.d.ethz.ch\share-images-1-$\');
        strNewPath = strrep(strNewPath,'/',filesep);
        strConvPath = strNewPath;
    elseif isunix
        strConvPath = strRootPath;
    end
    % and return
    return
end

% if no format is recognized on the input, throw warning and return (output
% is equal to input).
if ~intPathFormat
    
    % only throw warning if format_path was not in the stack trace
    ST = dbstack('-completenames');
    
    % end of function
    return
end

% determine available NAS paths on this machine for both share-2-$ and
% share-3-$
if ispc
    % add ".d.ethz.ch", if not present, strrep it out of here :-)
    strLocalBaseShare2 = '\\nas-unizh-imsb1.ethz.ch\share-2-$\';
    strLocalBaseShare3 = '\\nas-unizh-imsb1.ethz.ch\share-3-$\';
    strLocalBaseShare4 = '\\nas21nwg01.ethz.ch\biol_uzh_pelkmans_s4\';
    strLocalBaseShare5 = '\\nas21nwg01.ethz.ch\biol_uzh_pelkmans_s5\';
    strLocalBaseShare6 = '\\nas21nwg01.ethz.ch\biol_uzh_pelkmans_s6\';
    strLocalBaseShare7 = '\\nas21nwg01.ethz.ch\biol_uzh_pelkmans_s7\';
    strLocalBaseShare8 = '\\nas21nwg01.ethz.ch\biol_uzh_pelkmans_s8\';
elseif ismac
    if any(fileattrib('~/shares/'))
        try
            userDir = char(java.lang.System.getProperty('user.home')); % get user's home directory
        catch
            warning('Could not get user home directory. Check whether java is installed on your machine.\n')
            userDir = '';
        end
        strLocalBaseShare2 = sprintf('%s/shares/ethz-share2/',userDir);
        strLocalBaseShare3 = sprintf('%s/shares/ethz-share3/',userDir);
        strLocalBaseShare4 = sprintf('%s/shares/ethz-share4/',userDir);
        strLocalBaseShare5 = sprintf('%s/shares/ethz-share5/',userDir);
        strLocalBaseShare6 = sprintf('%s/shares/ethz-share6/',userDir);
        strLocalBaseShare7 = sprintf('%s/shares/ethz-share7/',userDir);
        strLocalBaseShare8 = sprintf('%s/shares/ethz-share8/',userDir);
    else
        strLocalBaseShare2 = '/Volumes/share-2-$/';
        strLocalBaseShare3 = '/Volumes/share-3-$/';
        strLocalBaseShare4 = '/Volumes/biol_uzh_pelkmans_s4/';
        strLocalBaseShare5 = '/Volumes/biol_uzh_pelkmans_s5/';
        strLocalBaseShare6 = '/Volumes/biol_uzh_pelkmans_s6/';
        strLocalBaseShare7 = '/Volumes/biol_uzh_pelkmans_s7/';
        strLocalBaseShare8 = '/Volumes/biol_uzh_pelkmans_s8/';
    end
elseif isunix
    strLocalBaseShare2 = '/BIOL/imsb/fs2/bio3/bio3/';
    strLocalBaseShare3 = '/BIOL/imsb/fs3/bio3/bio3/';
    strLocalBaseShare4 = '/BIOL/sonas/biol_uzh_pelkmans_s4/';
    strLocalBaseShare5 = '/BIOL/sonas/biol_uzh_pelkmans_s5/';
    strLocalBaseShare6 = '/BIOL/sonas/biol_uzh_pelkmans_s6/';
    strLocalBaseShare7 = '/BIOL/sonas/biol_uzh_pelkmans_s7/';
    strLocalBaseShare8 = '/BIOL/sonas/biol_uzh_pelkmans_s8/';
end

% some windowd machines need the ".d.ethz.ch", and some don't
if ispc   % note that on brutus: unmounted, but brutus is not pc
    if ~fileattrib(strLocalBaseShare2)
        strLocalBaseShare2 = '\\nas-biol-imsb-1.d.ethz.ch\share-2-$\';
        
        % check path with .d.ethz.ch, and resort to unizh address if not
        % present, as this is now the most likely available case
        if ~fileattrib(strLocalBaseShare2)
            strLocalBaseShare2 = strrep(strLocalBaseShare2,'.d.ethz.ch','');
        end
        
    end
    
    if ~fileattrib(strLocalBaseShare3)
        strLocalBaseShare2 = '\\nas-biol-imsb-1.d.ethz.ch\share-2-$\';
        
        % check path with .d.ethz.ch, and resort to unizh address if not
        % present, as this is now the most likely available case
        if ~fileattrib(strLocalBaseShare3)
            strLocalBaseShare3 = strrep(strLocalBaseShare3,'.d.ethz.ch','');
        end
    end
    
    
end

%%% check if they are both present, otherwise notify user (may get
%%% annoying) [TS: deactivated since on brutus share two not mounted:
%%% prevent excess fileattrib]
% if ~fileattrib(strLocalBaseShare2)
%     warning('naspathconv:unableToReachNasShare','%s: unable to reach share-2-$ on ''%s'', please map this path',mfilename,strLocalBaseShare2)
% end
% if ~fileattrib(strLocalBaseShare3)
%     warning('naspathconv:unableToReachNasShare','%s: unable to reach share-3-$ on ''%s'', please map this path',mfilename,strLocalBaseShare3)
% end
% Disabled warning for SoNas Shares (indeed got annoying)
% if ~fileattrib(strLocalBaseShare4)
%     warning('naspathconv:unableToReachSoNasShare','%s: unable to reach share-4 on ''%s'', please map this path',mfilename,strLocalBaseShare4)
% end
% if ~fileattrib(strLocalBaseShare6)
%     warning('naspathconv:unableToReachSoNasShare','%s: unable to reach share-6 on ''%s'', please map this path',mfilename,strLocalBaseShare6)
% end
% if ~fileattrib(strLocalBaseShare7)
%     warning('naspathconv:unableToReachSoNasShare','%s: unable to reach share-7 on ''%s'', please map this path',mfilename,strLocalBaseShare7)
% end


% determine the Nas Share number we want to get (0 is unknown/default)
intNasShare = 0;
strPathBase = char(regexpi(strRootPath,strShareRegexp,'Match'));
strNasShare = regexpi(strRootPath,strShareRegexp,'Tokens');
try
    intNasShare = str2double(strNasShare{1});
    if isempty(intNasShare)
        intNasShare = 0; % if not determinable, set to 0
    end
end
% if we weren't able to determine the share number, we might want to
% error-out (or alternatively try to determine the share number via
% trying...)
if (intNasShare == 0 | ~isnumeric(intNasShare)) & (intPathFormat > 1)
    error('%s: unable to determine share number from input path ''%s''',mfilename,strRootPath)
end



% do the conversion, by trying different nas shares
strNewPath = '';
if intNasShare == 0
    
    % create tentative path for share 2
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare2);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath2 = strrep(strNewPath,'\',filesep);
    
    % create tentative path for share 3
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare3);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath3 = strrep(strNewPath,'\',filesep);
    
    % create tentative path for share 4
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare4);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath4 = strrep(strNewPath,'\',filesep);
    
    % create tentative path for share 5
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare5);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath5 = strrep(strNewPath,'\',filesep);
    
    % create tentative path for share 6
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare6);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath6 = strrep(strNewPath,'\',filesep);
    
    % create tentative path for share 7
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare7);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath7 = strrep(strNewPath,'\',filesep);

    % create tentative path for share 8
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare8);
    strNewPath = strrep(strNewPath,'/',filesep);
    strNewPath8 = strrep(strNewPath,'\',filesep);

    
    % try SONAS shares at first, only then legacy share2/3, which are not
    % on brutus
    if fileattrib(strNewPath4)
        strNewPath = strNewPath4;
    elseif fileattrib(strNewPath6)
        strNewPath = strNewPath6;
    elseif fileattrib(strNewPath7)
        strNewPath = strNewPath7;
    elseif fileattrib(strNewPath8)
        strNewPath = strNewPath8;
    elseif fileattrib(strNewPath5)
        strNewPath = strNewPath5;
    elseif  fileattrib(strNewPath2)
        strNewPath = strNewPath2;
    elseif fileattrib(strNewPath3)
        strNewPath = strNewPath3;
    else
        % if we weren't able to determine the nas share, try to use ETH
        % Share 4
        warning('naspathconv:unableToDetermineNasShareNumber','%s: unable to determine the share-number from ''%s'', and neither share paths is found. Assuming share-4, but expext a crash..',mfilename,strPathBase)
        strNewPath = strNewPath4;
    end
    
    % do the conversion for share-2-$
elseif intNasShare == 2
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare2);
    % do the conversion for share-3-$
elseif intNasShare == 3
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare3);
    % do the conversion for share-4
elseif intNasShare == 4
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare4);
    % do the conversion for share-5
elseif intNasShare == 5
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare5);
    % do the conversion for share-6
elseif intNasShare == 6
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare6);
    % do the conversion for share-7
elseif intNasShare == 7
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare7);
elseif intNasShare == 8
    % do the conversion for share-8
    strNewPath = strrep(strRootPath,strPathBase,strLocalBaseShare8);
end
strNewPath = strrep(strNewPath,'/',filesep);
strNewPath = strrep(strNewPath,'\',filesep);

strConvPath = strNewPath;



end % function
% 
% function doNotWarnByNetwork = ignoreNetworkWarning
% % prevent that warning that path is not network path, will not be displayed
% % on some computers
% 
% 
% persistent cachedComputerToIgnore  % remove excess system calls by caching
% 
% if isempty(cachedComputerToIgnore)
%     
%     [~, name] = system('hostname');
%     name = lower(name);
%     
%     if strcmpi(name,'tswork'); % Thomas' computer
%         cachedComputerToIgnore = true;
%     elseif strcmpi(name,'pauli'); % Gabriele's computer
%         cachedComputerToIgnore = true;
%     else
%         cachedComputerToIgnore = false;
%     end
%     
% end
% 
% doNotWarnByNetwork = cachedComputerToIgnore;
% 
% end
% 
