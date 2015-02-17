% SHARED SETTINGS
S.OutputDirectory = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\Intensities');    % insert the name of the output folder
S.selChannel = 5;       % Channel of interest (note: you will likely have to adjust the function for metadata for your screening platform)
S.ZPositions = 1;       % Z Position of Interest, for consindering multiple layers use array, e.g.: 1:10 to include plane 1 to 10
S.fractionImages = 1;   % Fraction of images to be sampled
S.quantileOfMinimumIntensity = 0.01;      % Intensity quantile regarded as the minimum intensity of an image (note: use quantile and not minimum against unresponsive pixels)
S.quantileOfMaximumIntensity = 0.995;     % Intensity quantile regarded as the maximum intensity of an image (note: use quantile and not maximum against debris from plasticware)
S.applyIllumniationCorrection = true;     % apply NBBS illumination prior to analysis, options are true and false
S.downsampleFactor = 10;                  % Factor by which images are downsampled prior to determining intensities. Speeds up calcuation per image

% Define wells with positions of positive and negative controls. This
% follows standard regular expressions (which are used as a filter for the
% file names, that should be used for setting up analysis -> note this will
% depend upon the names from your microscope)
PosStandard = '_D18_|_insertNamesOfMultiplePosWells_'; % Hprt1 (a housekeeping gene expressed at intermediate levels)
NegStandard = '_D13_|_insertNamesOfMultipleNegWells_'; % dapB (a bacterial gene, that should be absent from human cells)

% Setup for each plate (note: here only one plate)
P{1}.name = 'MyExamplePlate';        
P{1}.pos = PosStandard;    
P{1}.neg = NegStandard;
P{1}.path = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExamplePlate\');   % insert the path to your plate of interest

% Settings for parallizing computation on different computational nodes in
% our cluster: You will likely have to adjust the commands and the
% "runDistributedJob.m" function for your cluster.
strBrutus = 'bsub -W 8:00'; % standard submission command for LSF based clusters
strFunction = 'SpotThrDetection.brutusIntensitiesOfPlate'; % function, which is called and actually obtains the intensities

% Create Batch file for tracing/parallelizing computational jobs on different nodes
numPlates = length(P);
for j=1:numPlates
    CurrName = P{j}.name;
    % add labels for human interpretabilty, note that in principle the script is not limited to positive and negative control wells
    P{j}.regImageGroups.description{1} = 'Hprt1 positive control wells'; %#ok<*SAGROW> 
    P{j}.regImageGroups.wells{1} = P{j}.pos;
    P{j}.regImageGroups.description{2} = 'dapB negative control wells';
    P{j}.regImageGroups.wells{2} = P{j}.neg;
end

SettingsFileName = fullfile(S.OutputDirectory,'Submission_SpotIntensities.mat');
ClusterSettings.Plate = P;
ClusterSettings.Shared = S;
ensurePresenceOfDirectory(S.OutputDirectory)
save(SettingsFileName,'ClusterSettings');

% Submit analysis of intensities into jobs parallelized by plate
for j=1:numPlates
    % note: this is submitted by a custom function, which can automatically
    % submit jobs to a cluster, if this code is executed on a cluster. You 
    % will likely have to adjust that function to support your specific cluster
    % environment  
    runDistributedJob(strBrutus,strFunction,SettingsFileName,j); 
end
