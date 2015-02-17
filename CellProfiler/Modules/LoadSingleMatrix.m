function handles = LoadSingleMatrix(handles)

% Help for the Load Single Matrix module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Loads a single matrix image, which will be used for all image cycles.
% *************************************************************************
% Note: for most purposes, you will probably want to use the Load Images
% module, not this one.
%
% similar to LoadSingeImage of orignal CellProfiler. 
% However, this module allows the import of relative image values higher than 1, 
% which can be useful for spot detection bias correction.
%
% See also LoadImages.

% 
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% *** Original Load Single Image Module ****
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
%
% Website: http://www.cellprofiler.org
%
% *** Load Single Matrix Module ****
% 
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%
% $Revision: 4533 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = This module loads one matrix as an image for *all* cycles that will be processed. Normally, use the Load Images module to load new sets of images during each cycle of processing.

%pathnametextVAR02 = Enter the path name to the folder where the images to be loaded are located.  Type period (.) for the default image folder.
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%filenametextVAR03 = What .mat file do you want to load? 
TextToFind{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call that image?
%defaultVAR04 = OrigBlue
%infotypeVAR04 = imagegroup indep
ImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% Remove slashes '/' from the input
tmp1 = {};
tmp2 = {};
for n = 1
    if ~strcmp(TextToFind{n}, 'NO FILE LOADED') && ~strcmp(ImageName{n}, 'Do not load')
        tmp1{end+1} = TextToFind{n};
        tmp2{end+1} = ImageName{n};
    end
end
TextToFind = tmp1;
ImageName = tmp2;

%%% Get the pathname and check that it exists
if strncmp(Pathname,'.',1)
    if length(Pathname) == 1
        Pathname = handles.Current.DefaultImageDirectory;
    else
        Pathname = fullfile(handles.Current.DefaultImageDirectory,Pathname(2:end));
    end
end
SpecifiedPathname = Pathname;
if ~exist(SpecifiedPathname,'dir')
    error(['Image processing was canceled in the ', ModuleName, ' module because the directory "',SpecifiedPathname,'" does not exist. Be sure that no spaces or unusual characters exist in your typed entry and that the pathname of the directory begins with / (for Mac/Unix) or \ (for PC).'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

if isempty(ImageName)
    error(['Image processing was canceled in the ', ModuleName, ' module because you have not chosen any images to load.'])
end

for n = 1:length(ImageName)
    %%% This try/catch will catch any problems in the load images module.
    try
        CurrentFileName = TextToFind{n};
        %%% The following runs every time through this module (i.e. for
        %%% every cycle).
        %%% Saves the original image file name to the handles
        %%% structure.  The field is named appropriately based on
        %%% the user's input, in the Pipeline substructure so that
        %%% this field will be deleted at the end of the analysis
        %%% batch.
        fieldname = ['Filename', ImageName{n}];
        handles.Pipeline.(fieldname) = CurrentFileName;
        fieldname = ['Pathname', ImageName{n}];
        handles.Pipeline.(fieldname) =  Pathname;

        FileAndPathname = fullfile(Pathname, CurrentFileName);
        %%%%% start %%%%%%%% [TS] %%%%%%%%%%%%%%%%      
        impImage = load(FileAndPathname, '-mat');
        impField = fieldnames(impImage);
        if length(impField) ==1
            LoadedImage = double(impImage.(impField{1}));
        end
              
        %[LoadedImage, handles] = CPimread(FileAndPathname,handles);
        
        %%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Saves the image to the handles structure.
        handles.Pipeline.(ImageName{n}) = LoadedImage;

    catch ErrorMessage = lasterr;
        ErrorNumber = {'first','second','third','fourth'};
        error(['Image processing was canceled in the ', ModuleName, ' module because an error occurred when trying to load the ', ErrorNumber{n}, ' set of images. Please check the settings. A common problem is that there are non-image files in the directory you are trying to analyze. Matlab says the problem is: ', ErrorMessage])
    end % Goes with: catch

    % Create a cell array with the filenames
    FileNames(n) = {CurrentFileName};
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The figure window display is unnecessary for this module, so the figure
%%% window is closed the first time through the module.
%%% Determines the figure number.
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
%%% Closes the window if it is open.
if any(findobj == ThisModuleFigureNumber)
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure('','NarrowText',ThisModuleFigureNumber)
    end
    for n = 1:length(ImageName)
        drawnow
        %%% Activates the appropriate figure window.
        currentfig=CPfigure(handles,'Text',ThisModuleFigureNumber);
        if iscell(ImageName)
            TextString = [ImageName{n},': ',FileNames{n}];
        else
            TextString = [ImageName,': ',FileNames];
        end
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 .85-(n-1)*.15 .95 .1],'BackgroundColor',[.7 .7 .9])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% First, fix feature names and the pathname
PathNames = cell(1,length(ImageName));
FileNamesText = cell(1,length(ImageName));
PathNamesText = cell(1,length(ImageName));
for n = 1:length(ImageName)
    PathNames{n} = Pathname;
    FileNamesText{n} = [ImageName{n}];
    PathNamesText{n} = [ImageName{n}];
end

%%% Since there may be several load modules in the pipeline which all write to the
%%% handles.Measurements.Image.FileName field, we have store filenames in an "appending" style.
%%% Here we check if any of the modules above the current module in the pipline has written to
%%% handles.Measurements.Image.Filenames. Then we should append the current filenames and path
%%% names to the already written ones.
if  isfield(handles,'Measurements') && isfield(handles.Measurements,'Image') && isfield(handles.Measurements.Image,'FileNames')
    if length(handles.Measurements.Image.FileNames) == SetBeingAnalyzed
        % Get existing file/path names. Returns a cell array of names
        ExistingFileNamesText = handles.Measurements.Image.FileNamesText;
        ExistingFileNames     = handles.Measurements.Image.FileNames{SetBeingAnalyzed};
        ExistingPathNamesText = handles.Measurements.Image.PathNamesText;
        ExistingPathNames     = handles.Measurements.Image.PathNames{SetBeingAnalyzed};

        % Append current file names to existing file names
        FileNamesText = cat(2,ExistingFileNamesText,FileNamesText);
        FileNames     = cat(2,ExistingFileNames,FileNames);
        PathNamesText = cat(2,ExistingPathNamesText,PathNamesText);
        PathNames     = cat(2,ExistingPathNames,PathNames);
    end
end

%%% Write to the handles.Measurements.Image structure
handles.Measurements.Image.FileNamesText                   = FileNamesText;
handles.Measurements.Image.FileNames(SetBeingAnalyzed)     = {FileNames};
handles.Measurements.Image.PathNamesText                   = PathNamesText;
handles.Measurements.Image.PathNames(SetBeingAnalyzed)     = {PathNames};