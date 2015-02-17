%%%%%%% General Settings  %%%%%%%%%%%

% indicate where intensities per image can be found (previously created by
% ObtainIntensities script)
strFolderWithIntensityMeasurements = nnpc(...
    '\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\Intensities');

strOutputFolder = nnpc(...
    '\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\Counts');

% References for finding absolute intensity thresholds based upon 
% Minima/Extrema of reference image sets will be set by quantiles of these
% measurments; We recommend to initially test different quantiles
Q = {[0.1 0.8 0.4 0.8];[0.8 0.8 0.6 0.8];}; % test multiple parameters for intensity limits

% This corresponds to the fraction of the distance between the upper 
% threshold for the minimal intensity and the lower for the maximal 
% intensity, which is used as an intensity threshold for the spots
% e.g. Intensity Thresholds are [100 200 400 1000] and sMin is 0.3 than 
% this would correspond to min pixel intensity of 260
MI = [0.05;0.3;];

useOldSettingsFile = false;

InvariantS.fractionImagesSpots = 1;                          % Fraction of images to be sampled
InvariantS.FilterMethod = '2D LoG';                            % Method used for generating the filter; see help of FSPECIALCP3D(). eg. '3D LoG, Raj' or '2D LoG'
InvariantS.FilterSize = 5;                                     % Size of Filter, approx. size of Object in Pixels.
InvariantS.SpotThresholdsToTest = 0.0001 : 0.002 :  0.3;       % Thresholds to be tested. Can be either a vector or a number.
InvariantS.ObjSizeThr = [];         % Minimal size of objects
InvariantS.CustomRescaleThr = [];   % Intensity Boundaries. See RESCALETHR in the help for OBJBYFILTER.
InvariantS.ObjIntensityThr = [];    % Miniamal intensity of pixels allowed to be part of a spot. See help of OBJBYFILETER()
InvariantS.closeHoles = true;       % note: closing holes becomes slow for 3D spot detection


useMinimalOfPositiveControlForMinThreshold = false; % optionally set to true, if plate is  dirty. in that case the lower intensity used for spot detection will be derived from positive control

% Settings for job submission
doOnlyTestRun = false; % if true, only one job should be sent out to the cluster (for testing if processing works)
strBrutus = 'bsub -W 8:00'; % standard submission command for LSF based clusters
strFunction = 'SpotThrDetection.brutusSpotsOfPlateSeparateIllCorr'; % function which obtains spot counts per image

% Compute additional shared settings
InvariantS.Filter = fspecialCP3D(InvariantS.FilterMethod,InvariantS.FilterSize);
InvariantS.InputDirectory = strFolderWithIntensityMeasurements;
InvariantS.InputFileBase = 'SpotSetupIntensities'; % default name of output of previous script

BaseOfIntensitySubmission  = 'Submission_SpotIntensities';
nameOfIntensitySettingsFile = fullfile(InvariantS.InputDirectory,getFilesAndDirectories(InvariantS.InputDirectory, BaseOfIntensitySubmission));
if numel(nameOfIntensitySettingsFile) == 1
    IntensitySettingsFile = loadd(nameOfIntensitySettingsFile{:});
else
    error('potentially no or multiple intensity settings files')
end

numPlatesToAnalyse = length(IntensitySettingsFile.Plate);   % get information about plates to process
pathToImportFile = cell(1,numPlatesToAnalyse);
for j=1:numPlatesToAnalyse
   P{j}.name = IntensitySettingsFile.Plate{j}.name;  %#ok<*SAGROW>  % for traceability of analysis
   pathToImportFile{j} = fullfile(InvariantS.InputDirectory,[InvariantS.InputFileBase P{j}.name '.mat']);
end

%%%%% Computation of spots  %%%%%
ensurePresenceOfDirectory(strOutputFolder);
doRun = true; % internal variable, whether to continue (needed for testrun option)

if useMinimalOfPositiveControlForMinThreshold == true
    STRINGuseMinimalOfPositiveControlForMinThreshold = 'true';
elseif useMinimalOfPositiveControlForMinThreshold == false
    STRINGuseMinimalOfPositiveControlForMinThreshold = 'false';
else
    error('did not recognize format of useMinimalOfPositiveControlForMinThreshold, must be either true or false');
end

for jq = 1:length(Q) % loop through testing different intensity quantiles
    for jmi = 1: length(MI) % loop through testing different minimal intensities
        
        % build settings file for one set of spot detection settings, that should
        % be tested. Later parallelization will be done by single plates
        S = InvariantS; % Use invariant settings as a template
        S.vSelectionIntensityQuantiles = Q{jq};
        S.sMinIntensity = MI(jmi);
        
        % make name of output sub-folder, which contains information on varied
        % submission settings
        CurrSubfolderName = ['SpotCount' num2str(S.vSelectionIntensityQuantiles) 'l' num2str(MI(jmi))];
        CurrSubfolderName = regexprep(CurrSubfolderName, '         ', 'd');
        posPoint = strfind(CurrSubfolderName,'.');
        CurrSubfolderName(posPoint) = 'p';
        
        % make output sub-folder
        S.OutputDirectory = nnpc(fullfile(strOutputFolder,CurrSubfolderName));    % folder where output should saved
        ensurePresenceOfDirectory(S.OutputDirectory);
        
        S.OutputName = ['SpotSetupCount' CurrSubfolderName];
               
        %%%%%%%%%%%%%%%	JOB SUBMISSION %%%%%%%%%%%%%%%%%%%%%%%%%
        % Save Settings for spot detections to Separate File
        
        S.SubmissionFile = fullfile(S.OutputDirectory,'Submission_SpotCount.mat');
        
        ClusterSettings.Plate = P;
        ClusterSettings.Shared = S;
        
        if any(fileattrib(S.SubmissionFile));
            fprintf('Do not submit, as jobs submitted once \n');    % here conservative; might change depending on your policy / job scheduling framework
            continue;
        else
            save(S.SubmissionFile,'ClusterSettings');
        end
        
        % Parallelize by plate for each variation of the spot detection
        % settings
        for PlateIndex=1:numPlatesToAnalyse
            if doRun == true
                % note: this is submitted by a custom function, which adjusts jobs
                % according to our different compuational requirements
                runDistributedJob(strBrutus,strFunction,S.SubmissionFile,pathToImportFile{PlateIndex}, PlateIndex, STRINGuseMinimalOfPositiveControlForMinThreshold);
            end
            
            if doOnlyTestRun == true
                doRun = false;
            end
            
        end
        
    end
end