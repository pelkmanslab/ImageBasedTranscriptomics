%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%

% Data to import
CurrentDetectionVersion = 'SpotCount0p8d0p8d0p6d0p8l0p3';  % name of subfolder which represennts the quantiles used for reference bounds of intensities
pathImport = nnpc(fullfile('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\Counts',...
    CurrentDetectionVersion));

% Visualization 
ZoomOnYSpotCount = [20000 15000 1000 50]; % Zoom for Absolute Spot Number. Options: [] or number or vector of arbitrary length; Subfigures are created for all provided numbers.
ZoomOnYRelativeCount = [2 1.2 1];


%%%% CALCULATION %%%%%%

baseImport = ['SpotSetupCount' CurrentDetectionVersion];

% get Identifies of plates 
strBatchFile = fullfile(pathImport, 'Submission_SpotCount.mat');
BatchFile = loadd(strBatchFile);
numPlates = length(BatchFile.Plate);
for j=1:numPlates;
    P{j}.name = BatchFile.Plate{j}.name;
end

% Read results
fileExists = false(1,length(P));
for k=1:numPlates
    pathToFile = fullfile(pathImport,[baseImport P{k}.name '.mat']);
    try
        load(pathToFile);
        fileExists(k) = true;
    catch CouldnotLoadFile
        fprintf([P{k}.name ' was not found.\n']);
    end
    rescalingQuantiles{k} = [   strSpotCount.PriorInputSettings.Shared.quantileOfMinimumIntensity, ...
        strSpotCount.PriorInputSettings.Shared.quantileOfMaximumIntensity]; %#ok<*SAGROW>
    ObjCount{k} = strSpotCount.Output.ObjCount;
    ObjIntensityThr{k} = strSpotCount.Output.ObjIntensityThr;
    vIntensityBoundaries{k} = strSpotCount.Output.vIntensityBoundaries;
    SpotThresholdsToTest{k} = strSpotCount.InputSettings.Shared.SpotThresholdsToTest;
    FilterSize{k} = strSpotCount.InputSettings.Shared.FilterSize;
end

%%%% CALCULATION %%%%%%
fprintf('Please manually choose threshold, see ValueChosenThreshold in script \n');

k=1; % plate of interest
ValueChosenThreshold{k}=0.07;  % please insert threshold, which you consider good);
SpotThrDetection.visualizeSpotThresholds(ObjCount{k}, SpotThresholdsToTest{k}, 'Absolute', ValueChosenThreshold{k},ZoomOnYSpotCount); title(P{k}.name);
SpotThrDetection.visualizeSpotThresholds(ObjCount{k}, SpotThresholdsToTest{k}, 'Relative', ValueChosenThreshold{k},ZoomOnYRelativeCount);title(P{k}.name);

% To see all the other settings used for spot detection, have a look at
% following variables (or - not shown: - programmaticaly include them in a 
% CellProfiler pipeline)
FilterSize{k}
rescalingQuantiles{k}
vIntensityBoundaries{k}
ObjIntensityThr{k} 

