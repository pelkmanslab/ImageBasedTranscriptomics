function [matMeanImage matStdImage hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,cacheInRam)
if nargin < 3
    cacheInRam = false;
end

strPathToCurrentIllumination = fullfile(strBatchDir,...
    sprintf('Measurements_batch_illcor_channel%03d_zstack000.mat',iChannel));

matMeanImage = [];
matStdImage =[];

if ~any(fileattrib(strPathToCurrentIllumination))
    hasIlluminationCorrection = false;
    warning('matlab:bsBla','%s:  failed to load illumination correction %s',mfilename,strPathToCurrentIllumination);
else
    hasIlluminationCorrection = true;
    if cacheInRam == false
        ImportedIlluminationCorrection = load(strPathToCurrentIllumination);
        matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
        matStdImage = double(ImportedIlluminationCorrection.stat_values.std);
    else
        ImportedIlluminationCorrection = cacheForIlluminationStats(strPathToCurrentIllumination);
        matMeanImage = ImportedIlluminationCorrection.stat_values.mean;  % conversion to double in cache to save time
        matStdImage = ImportedIlluminationCorrection.stat_values.std;
    end
    
end

end

function dat = cacheForIlluminationStats(strFileName)
% Initialize Persistent variables for caching
persistent CachedMeasurments;
persistent OriginalPathOfChached;

if isempty(CachedMeasurments)
    CachedMeasurments = cell(0);
end

if isempty(OriginalPathOfChached)
    OriginalPathOfChached = cell(0);
end

nStrFileName = npc(strFileName); % npc to ensure that each file only stored once in cache once

[isCached cachedLoc]= ismember(nStrFileName,OriginalPathOfChached);

if ~isCached   % load into cache, if absent there
    fprintf('Caching illumination correction ... ');
    cachedLoc = length(CachedMeasurments) + 1;
    
    ImportedIlluminationCorrection = load(nStrFileName);
    ImportedIlluminationCorrection.stat_values.matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
    ImportedIlluminationCorrection.stat_values.std = double(ImportedIlluminationCorrection.stat_values.std);
    
    OriginalPathOfChached{cachedLoc} = nStrFileName;
    CachedMeasurments{cachedLoc} = ImportedIlluminationCorrection;    
    
    fprintf('complete \n');
end

dat = CachedMeasurments{cachedLoc}; % retreive data
end