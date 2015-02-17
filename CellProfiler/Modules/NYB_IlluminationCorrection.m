function handles = NYB_IlluminationCorrection(handles)
% Help for the "Nothing Yields Better" Ilumination Correction module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates a mean illumination function and the std of the illumination
% mean, used to correct uneven illumination/lighting/shading via Z-scoring.
% *************************************************************************
%
% This module calculates the mean and std values for the intensity of each
% pixel position in a set of images from the same wave length. Such
% fuctions can later be used to to correct the intensity values of each
% pixel by substracting the mean function and normalizing with the std
% function. Such proceedure sould not only normalize for the background
% intensities but also for amplitude differences whic result from optical
% aberrances or differences in illumination.
% Note: this module uses only a precalculated illumination correction
% files (from iBrain) that can be placed in the default output folder.
%
% Authors:
%   Nico Battich
%   Berend Snijder
%   Yauhen Yakimovich
%
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images to be used to calculate the illumination functions (input)?
%infotypeVAR01 = imagegroup
InputName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How do you want to call the corrected images (output)?
%defaultVAR02 = CorrBlue
%infotypeVAR02 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,2});


%textVAR03 = Do you also want to smooth mean/std statistics values before using them for correction?
%choiceVAR03 = No
%choiceVAR03 = Yes
DoMeanStdSmoothing = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Smoothing is done by a Gaussian filter to both illumination functions. Enter the size of the filter.
%defaultVAR04 = 100
SmoothingSize = str2num(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Enter the factor 'x' for the sigma calculation of the Gaussian filter. where sigma = x*size.
%defaultVAR05 = 0.5
SmoothingSigma = str2num(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Enter the relative path name to the folder where the illumination correction files are located (starting with "./../"). Type "Pre" (Pre) to load files from previous multiplexing cycle. Type period (.) for default directory.
%defaultVAR06 = .
AlternativeIllCorrFolder = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%%%VariableRevisionNumber = 4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% The illumination correction image was calculated using all the incoming
% images by iBRAIN. Load the mean and std statistics assuming iBrain file 
% structure. 


% Get the channel of the image.
[intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});
intZstackNumber = 0;
disp(sprintf('Applying illumination correction on channel number: %d', intChannelNumber))

% Store the stat in the Measurement field, with channel and szstack
% specific fieldnames.
if strcmp(AlternativeIllCorrFolder,'Pre')
   strStatFieldName = sprintf('illcor_ch%03dz%03d_Pre',intChannelNumber,intZstackNumber);
else
   strStatFieldName = sprintf('illcor_ch%03dz%03d',intChannelNumber,intZstackNumber);
end

if handles.Current.SetBeingAnalyzed == 1
    % Get the BATCH directory and load illcor statistics.
    
    % Get path to folder where files are located
    if strncmp(AlternativeIllCorrFolder,'.',1)
        if length(AlternativeIllCorrFolder) == 1
            VarPathname = handles.Current.DefaultOutputDirectory;
        else
            VarPathname = fullfile(handles.Current.DefaultOutputDirectory,AlternativeIllCorrFolder(2:end));
            fprintf(['=================================================' ...
                     '=================================================' ...
                     '=================================================' ...
                     '\n\n%s: You specified an alternative path:\n%s\n\n' ...
                     '=================================================' ...
                     '=================================================' ...
                     '=================================================' ...
                     '\n'],mfilename,VarPathname);        
        end
    elseif strcmp(AlternativeIllCorrFolder,'Pre')
        if isfield(handles.Pipeline,['Pathname',InputName])
            VarPathname = fullfile(strrep(handles.Pipeline.(['Pathname',InputName]),'TIFF','BATCH'));
            fprintf(['=================================================' ...
                     '=================================================' ...
                     '=================================================' ...
                     '\n\n%s: You specified an alternative path in the LoadImages module:\n%s\n\n' ...
                     '=================================================' ...
                     '=================================================' ...
                     '=================================================' ...
                     '\n'],mfilename,VarPathname);
        else
            error(['Image processing was canceled in the ', ModuleName, ' module because the fieldname "',['Pathname',InputName],'" does not exist within handles.Pipeline. Be sure that you have saved it correctly using the LoadImages module'])
        end
    end
    strBatchDir = VarPathname;
    if ~exist(strBatchDir,'dir')
        error(['Image processing was canceled in the ', ModuleName, ' module because the directory "',strBatchDir,'" does not exist. Be sure that no spaces or unusual characters exist in your typed entry and that the pathname of the directory begins with /.'])
    end
    TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));    
    
    TempStats.stat_values.mean = double(TempStats.stat_values.mean);
    TempStats.stat_values.std = double(TempStats.stat_values.std);

    % Optional smoothing of illumination statistics computed by iBRAIN.
    switch DoMeanStdSmoothing
        case 'Yes'
            % Create gaussian filter handle for for smoothing
            H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);
            handles.Measurements.Image.([strStatFieldName,'_mean']) =  imfilter(TempStats.stat_values.mean,H,'replicate');
            handles.Measurements.Image.([strStatFieldName,'_std']) = imfilter(TempStats.stat_values.std,H,'replicate');
        case 'No'
            handles.Measurements.Image.([strStatFieldName,'_mean']) = TempStats.stat_values.mean;
            handles.Measurements.Image.([strStatFieldName,'_std']) = TempStats.stat_values.std;
    end
    clear TempStats;
end


%%%%%%%%%%%%%%%%%%
%%% LOAD IMAGE %%%
%%%%%%%%%%%%%%%%%%
strImageToImport = fullfile( ...
    handles.Pipeline.(strcat('Pathname',InputName)), ...
    handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});
OrigImage = double(imread(strImageToImport));
         

%%%%%%%%%%%%%%%%%%
%%% CORRECTION %%%
%%%%%%%%%%%%%%%%%%
IllumFilt_Mean = handles.Measurements.Image.([strStatFieldName,'_mean']);
IllumFilt_STD = handles.Measurements.Image.([strStatFieldName,'_std']);

% Avoid -Inf values after log10 transform.
OrigImage(OrigImage == 0) = 1;
% Apply z-score normalization for each single pixel.
ImageOutput = (log10(OrigImage) - IllumFilt_Mean) ./ IllumFilt_STD;
% Reverse z-score.
ImageOutput = (ImageOutput .* mean(IllumFilt_STD(:))) + mean(IllumFilt_Mean(:));
% Reverse log10 transform that was applied to images when learning 
% mean/std statistics as well the corrected image.
ImageOutput = 10 .^ ImageOutput;

% store non-scaled for visualization
ImageOutputPlot = ImageOutput;
% Rescale from 0  to 1
ImageOutput = ImageOutput/65535;
% fix potentially bad pixels (eg. dead pixels: std of 0)
ImageOutput = fixNonNumericalValueInImage(ImageOutput);
%save to handle structure
handles.Pipeline.(OutputName) = ImageOutput;
fieldname = ['Filename', OutputName];
handles.Pipeline.(fieldname) = {handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed}};
fieldname = ['Pathname', OutputName];
handles.Pipeline.(fieldname) = handles.Pipeline.(strcat('Pathname',InputName));


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

if any(findobj == ThisModuleFigureNumber)
       
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    subplot(2,2,1);    
   
    CPimagesc(IllumFilt_Mean,handles);
    colormap('JET')
    colorbar
    title('Mean Intensity Filter [Log10(intensity)]')
    
    subplot(2,2,3);    
    CPimagesc(IllumFilt_STD,handles);
    colormap('JET')
    colorbar
    title('STD Intensity Filter [Log10(intensity)]')    

    subplot(2,4,3);    
    CPimagesc(OrigImage,handles);
    colormap('JET')
    title('Original Image')
    
    subplot(2,4,4);
    CPimagesc(ImageOutputPlot,handles);
    colormap('JET')
    title('Corrected Image')    

    subplot(2,4,7);    
    hold on
    hist(OrigImage(:),round(length(OrigImage(:))/20))    
    hold off
    set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
    ylabel('Pixel Count')
    xlabel('Intensity')
    title('Original Image Histogram')

    subplot(2,4,8);
    hist(ImageOutputPlot(:),round(length(ImageOutputPlot(:))/20))
    set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
    title('Corrected Image Histogram')    
    ylabel('Pixel Count')
    xlabel('Intensity')
    
    drawnow
end    
