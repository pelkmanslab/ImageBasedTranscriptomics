function Images = imreadCP3D(strFilenames,datatype,illum_stat_values)
%IMREADCP3D generates a matrix containing the images specified in the
%   the array STRFILENAMES, where index positions in STRFILENAMES
%   correspond to different Z-planes. STRFILENAMES has to contain full
%   paths to the file to load.
%
%   Optionally the datatype of the imported image information can be
%   selected as 'double' 'single' or 'uint16'. If datatype is not specified
%   or [], 'double' will be used as default. Note that this option allows
%   to reduce the amount of memory required during the loading of multiple
%   images as well as the memory for their storage.
%   ----------------------------------
%   Information about module:
%   Created for CP3D to load multiple images into one stack
%
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHECK INPUT  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if datatype has been selected
if nargin < 2
    selDataType = 'double';
elseif isempty(datatype)
    selDataType = 'double';
else
    selDataType = (datatype);
end

if nargin < 3
    doIllumCorrection = false;
else
    doIllumCorrection = true;
    checkInputOfStatValues(illum_stat_values);
end

% Check if strFilenames are unambiguous
if size(strFilenames,1)>size(strFilenames,2)
    strFilenames = strFilenames';
end

if size(strFilenames,1) > 1
    error('strFilenames must be array of Filenames with n x 1 dimensions');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFiles = size(strFilenames,2);

% Initialize Variable for output image
% read image header of first image of one series to obtain rows and column
% dimensions
info = imfinfo(strFilenames{1,1});
Dim(1) = info.Height;
Dim(2) = info.Width;
Dim(3) = numFiles;

switch selDataType    % initialize according to datatype.
    case 'double'
        Images = zeros(Dim,'double');
    case 'single'
        Images = zeros(Dim,'single');
    case 'uint16'
        Images = zeros(Dim,'uint16');
    otherwise
        error('Datatype for loading images is not supported')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  LOAD IMAGES INTO MATRIX %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:numFiles
    try
        ImportedRawImage = double(imread(strFilenames{1,k})); % import as double for illumination correction
    catch CanNotLoad
        error(['Could not load image ' strFilenames{1,k} ' . Please check if file exists or if file format is supported (such as .tif or .png).']);
    end
    
    if doIllumCorrection == true
        isLog = 1; % [TS 140808] List here for clarity if output of iBrain statistics change in future.
        matMean = double(illum_stat_values.mean);
        matStd = double(illum_stat_values.std);
        ImageForLayer = IllumCorrect(ImportedRawImage,matMean,matStd,isLog);
    else
        ImageForLayer = ImportedRawImage;
    end
    
    switch selDataType
        case 'double'
            Images(:,:,k) = ImageForLayer;
        case 'single'
            Images(:,:,k) = ImageForLayer;
        case 'uint16'
            Images(:,:,k) = ImageForLayer;
        otherwise
            error('Datatype for loading images is not supported')
    end
    
    
end



end

function checkInputOfStatValues(illum_stat_values)
if ~isfield(illum_stat_values,'mean')
    error('Illumination correction file does not have mean defined');
elseif ~isfield(illum_stat_values,'std')
    error('Illumination correction file does not have mean defined');
end
end