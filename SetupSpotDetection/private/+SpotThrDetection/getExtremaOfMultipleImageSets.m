function [minIntensity, maxIntensity, ImageNames,IllumCorrRef] = getExtremaOfMultipleImageSets(PlatePath,regImageSets,selChannel,ZPositions,fractionImages,quantileMin,quantileMax,bnApplyIlluminationCorrection,downsampleFactor)
% Will create cells MININTENSITY and MAXINTENSITY, which lists for
% FRACTIONIMAGS randomly selected images of each
% REGIMAGESETS/ZPOSITIONS/SELCHANNEL the minimal
% and maximal intensity as specified by QUANTILEMIN and QUANTILEMAX.
% Quantiles are used to prevent that object detection is strongly
% influenced by outliers (note that CV7k has one false-active pixel with 10000).
%
% BNAPPLYILLUMINATIONCORRECTION can be true or false and defines, whether
% illumination correction should be performed
%
% DOWNSAMPLEFACTOR describes a factor by which images should be downsampled
% before indentifying Minima and Maxima. If it is too small, the
% caluclation of the Minima and Maxima will become the time-limiting factor
% of finding Thresholds for spots. Should be at least 10 for Yokogawa
% images (10 Megapixel)
%
% The outputs IMAGENAMES and ILLUMCORRREF are also picked up
% by sibling functions, which work together to finding spot references.
% They are outputted to prevent recalculation/loading to save time and
% prevent inconsistencies.
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html

% Initialize Image Sets according to the specified filters
numImageGroups = length(regImageSets);
ImageNames = cell(1,numImageGroups);
try
    for k = 1:numImageGroups
        ImageNames{k} = getFilteredMediafilenamesCP3D(fullfile(PlatePath,'/TIFF'),regImageSets{k},selChannel,ZPositions);
    end
catch CouldNotFindImageFiles
    error('Could not find files matching the description. Check name of channel and wells of interest')
end

% If requested, Load Illumination correction. Note that Ill correction
% reference images will be kept in memory to prevent reloading; Note that
% the following lines assume an equal illumination function for all z
if bnApplyIlluminationCorrection == true
    % Note: the following lines are the ones of original code used in
    % analysis of Hela; They are shown for reproduciblity. They were
    % replaced by new code which accounts for the change in data format of
    % illumination statistcs
    
    %     TempStats = load(fullfile(PlatePath,'/BATCH',sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',selChannel,0)));
    %     IllumFilt_Mean = TempStats.stat_values.mean;
    %     IllumFilt_STD = TempStats.stat_values.std;
    %     clear TempStats;   % clear since Illumination correction file contains information, which is not required / save memory
    
    strBatchDir = fullfile(PlatePath,'/BATCH');
    
    [IllumFilt_Mean IllumFilt_STD] = getIlluminationReference(strBatchDir,selChannel);
    
else
    IllumFilt_Mean = [];
    IllumFilt_STD = [];
end

PresentImageSets = cell2mat(cellfun(@(x) size(x,1), ImageNames, 'UniformOutput', false));
NumImagesToAnalye = ceil(PresentImageSets.*fractionImages);

minIntensity = cell(1,numImageGroups);
maxIntensity = cell(1,numImageGroups);

% Find Minimal and Maximal Intensity of images (as specified by Quantile)
for l=1:numImageGroups
    randIX = randperm(size(ImageNames{l},1));
    minIntensity{l} = nan(NumImagesToAnalye(l),1);
    maxIntensity{l} = nan(NumImagesToAnalye(l),1);
    
    for k = 1:NumImagesToAnalye(l);
        [minIntensity{l}(k) maxIntensity{l}(k)] = ...
            getImageIntensityExtremaCP3D(fullfile(PlatePath,'/TIFF/',ImageNames{l}{randIX(k)}),...
            quantileMin,quantileMax,downsampleFactor,IllumFilt_Mean,IllumFilt_STD);
        
        
        FractionAnalyzsed = k./NumImagesToAnalye(l).*100;
        if mod(k,floor(NumImagesToAnalye(l)./10))==0
            fprintf('%0.0f%% of image set %d of %d \n',FractionAnalyzsed,l,numImageGroups);
        end
    end
end

% Safe Reference of illumination correction to prevent reloading for
% downstream functions / ensure that correction is done by same reference
IllumCorrRef.IllumFilt_Mean = IllumFilt_Mean;
IllumCorrRef.IllumFilt_STD = IllumFilt_STD;

end