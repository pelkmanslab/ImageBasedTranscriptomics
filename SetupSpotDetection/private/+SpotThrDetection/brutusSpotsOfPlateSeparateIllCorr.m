function brutusSpotsOfPlateSeparateIllCorr(SettingsFile,IntensityFile,PlateIndex,useMinimalOfPositiveControlForMinThreshold,saveIlluminationCorrection)

if ischar(PlateIndex)
    PlateIndex = str2num(PlateIndex);
end

if nargin < 4
    useMinimalOfPositiveControlForMinThreshold = false;
else
    if ischar(useMinimalOfPositiveControlForMinThreshold)
        useMinimalOfPositiveControlForMinThreshold = str2num(useMinimalOfPositiveControlForMinThreshold); %#ok<ST2NM>
    end
end



if nargin < 5
    saveIlluminationCorrection = false;
else
    if ischar(saveIlluminationCorrection)
        saveIlluminationCorrection = str2num(saveIlluminationCorrection); %#ok<ST2NM>
    end
end



load(nnpc(SettingsFile));
load(nnpc(IntensityFile));

% Determine which settings should be used
prevID = handles.CurrentPlateIndex;
OutputDirectory = ClusterSettings.Shared.OutputDirectory;
OutputName = ClusterSettings.Shared.OutputName;

PlateName = ClusterSettings.Plate{PlateIndex}.name
if PlateName ~= handles.InputSettings.Plates{prevID}.name
    error(['Names of Plates do not match' PlateName handles.InputSettings.Plates{prevID}.name]);
end
selSpecifiedPathname = handles.InputSettings.Plates{prevID}.path;
quantileOfMinimumIntensity = handles.InputSettings.Shared.quantileOfMinimumIntensity;
quantileOfMaximumIntensity = handles.InputSettings.Shared.quantileOfMaximumIntensity;

% read output of previous function (which was determining intensities)
IllumCorrRef = handles.Output.IllumCorrRef;
ImageNames = handles.Output.ImageNames;
minIntensity = handles.Output.minIntensity;
maxIntensity = handles.Output.maxIntensity;

% load shared parameters for spot detection
Filter = ClusterSettings.Shared.Filter;
vSelectionIntensityQuantiles =  ClusterSettings.Shared.vSelectionIntensityQuantiles;
fractionImagesSpots =           ClusterSettings.Shared.fractionImagesSpots;
sMinIntensity =                 ClusterSettings.Shared.sMinIntensity;
CustomRescaleThr =              ClusterSettings.Shared.CustomRescaleThr;
ObjIntensityThr =               ClusterSettings.Shared.ObjIntensityThr;
SpotThresholdsToTest =          ClusterSettings.Shared.SpotThresholdsToTest;
closeHoles  =                   ClusterSettings.Shared.closeHoles;
ObjSizeThr =                    ClusterSettings.Shared.ObjSizeThr;


joinedName = [OutputName PlateName];
outputFilename = fullfile(...
    OutputDirectory,...
    [joinedName '.mat']);

% Obtain Intensity Thrshold(s)
% - for rescaling images
if isempty(CustomRescaleThr)
    if useMinimalOfPositiveControlForMinThreshold == false % default
        
        [vIntensityBoundaries] = SpotThrDetection.suggestImageIntensityBoundaries(...
            minIntensity{2},vSelectionIntensityQuantiles(1),...
            maxIntensity{2},vSelectionIntensityQuantiles(2),...
            maxIntensity{1},vSelectionIntensityQuantiles(3),...
            maxIntensity{1},vSelectionIntensityQuantiles(4));
        
    else   %  use minimal intensities of positive control as reference for spot detection
        [vIntensityBoundaries] = SpotThrDetection.suggestImageIntensityBoundaries(...
            minIntensity{1},vSelectionIntensityQuantiles(1),...
            minIntensity{1},vSelectionIntensityQuantiles(2),...
            maxIntensity{1},vSelectionIntensityQuantiles(3),...
            maxIntensity{1},vSelectionIntensityQuantiles(4));
        
    end
else
    vIntensityBoundaries = CustomRescaleThr;
    fprintf('Using custom intensity thresholds instead of recommended ones.');
end
% - minimal intensity of pixels
if isempty(ObjIntensityThr)
    if ~isempty(sMinIntensity)
        ObjIntensityThr = vIntensityBoundaries(2) + sMinIntensity*(vIntensityBoundaries(3)-vIntensityBoundaries(2));
    else
        fprintf('Using custom minimal intensity for pixels within spot instead of recommended ones.');
    end
end


%%% Save the graphical representation of the chosen threshold %%%

%SpotThrDetection.visualizeIntensityThresholds(minIntensity, maxIntensity, vIntensityBoundaries,[80 240])
%myTitle = [joinedName num2str(vSelectionIntensityQuantiles)];
%title(myTitle);
%gcf2pdf(OutputDirectory,myTitle);


%%%% Do the spot detection %%%%%

[ObjCount] = SpotThrDetection.runObjByFilterOnMultipleImages(...
    ImageNames,selSpecifiedPathname,Filter,SpotThresholdsToTest,[quantileOfMinimumIntensity quantileOfMaximumIntensity],vIntensityBoundaries,...
    ObjIntensityThr,closeHoles,ObjSizeThr,IllumCorrRef,fractionImagesSpots);


%%%%% SAVE OUTPUT %%%%%%%%%


strSpotCount.PriorPlateIndex = handles.CurrentPlateIndex;
strSpotCount.PriorInputSettings = handles.InputSettings;
strSpotCount.Input.minIntensity = minIntensity;
strSpotCount.Input.maxIntensity = maxIntensity;
strSpotCount.Input.ImageNames = ImageNames;

strSpotCount.CurrentPlateIndex = PlateIndex;
strSpotCount.InputSettings.Shared = ClusterSettings.Shared;
strSpotCount.InputSettings.Plates = ClusterSettings.Plate;
strSpotCount.Output.ObjCount = ObjCount;
strSpotCount.Output.ObjIntensityThr = ObjIntensityThr;
strSpotCount.Output.vIntensityBoundaries = vIntensityBoundaries; %#ok<STRNU>


save(outputFilename,'strSpotCount')

outputFilename = fullfile(...
    OutputDirectory,...
    [joinedName '_Illumination.mat']);

if saveIlluminationCorrection == true
    strSpotCountSupport.Input.IllumCorrRef = IllumCorrRef; %#ok<STRNU>
    save(outputFilename,'strSpotCountSupport')
end

fprintf('%s: plate %s completed successfully \n', mfilename, PlateName);


end