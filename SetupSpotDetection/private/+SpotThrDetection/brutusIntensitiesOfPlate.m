function brutusIntensitiesOfPlate(SettingsFile,PlateIndex)
load(SettingsFile);

if ischar(PlateIndex)
    PlateIndex = str2num(PlateIndex);
end

% Determine which settings should be used
OutputDirectory = ClusterSettings.Shared.OutputDirectory;
PlateName = ClusterSettings.Plate{PlateIndex}.name
selSpecifiedPathname = ClusterSettings.Plate{PlateIndex}.path;
regImageGroups = ClusterSettings.Plate{PlateIndex}.regImageGroups.wells;
selChannel = ClusterSettings.Shared.selChannel;
ZPositions = ClusterSettings.Shared.ZPositions;
fractionImagesIntensities = ClusterSettings.Shared.fractionImages;
quantileOfMinimumIntensity = ClusterSettings.Shared.quantileOfMinimumIntensity;
quantileOfMaximumIntensity = ClusterSettings.Shared.quantileOfMaximumIntensity;
applyIllumniationCorrection = ClusterSettings.Shared.applyIllumniationCorrection;
downsampleFactor = ClusterSettings.Shared.downsampleFactor;

% Obtain Intensity measurments
[minIntensity, maxIntensity, ImageNames, IllumCorrRef] = ...
    SpotThrDetection.getExtremaOfMultipleImageSets(...
    selSpecifiedPathname,regImageGroups,selChannel,ZPositions,fractionImagesIntensities,...
    quantileOfMinimumIntensity,quantileOfMaximumIntensity,...
    applyIllumniationCorrection,downsampleFactor);

% Save output together with all settings (for reproducibility)
handles.CurrentPlateIndex = PlateIndex;
handles.InputSettings.Shared = ClusterSettings.Shared;
handles.InputSettings.Plates = ClusterSettings.Plate;
handles.Output.minIntensity = minIntensity;
handles.Output.maxIntensity = maxIntensity;
handles.Output.ImageNames = ImageNames;
handles.Output.IllumCorrRef = IllumCorrRef;

outputFilename = fullfile(...
    OutputDirectory,...
    ['SpotSetupIntensities' ...
    PlateName '.mat']);

save(outputFilename,'handles')

fprintf('%s: job %s completed successfully', mfilename, ClusterSettings.Plate{PlateIndex}.name);
end