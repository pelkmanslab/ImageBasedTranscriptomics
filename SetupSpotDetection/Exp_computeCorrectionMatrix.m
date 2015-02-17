fprintf('Please note that, once the correction matrix is computed, you have to move the Correction Matrix to the BATCH folder of the plate, where it should be applied \n');

%%%%%% INPUT PARAMETERS  %%%%%%%%%%%%%%%
S.SpotSmoothing = 150;          % size of smoothing the coordinates of spots (pixel)
S.CorrectionSmoothing = 150;    % size of smoothing the correction matrix (pixel)
S.ImageDimensions = [2160 2560];    % [rows columns] of single images obtained by your microscope (pixel)
S.scaleToInputThreshold = true;     % scale correction matrix such that it uses the previously chosen references threshold (which must be input in CP module ScanSpotThresholds.m)

S.PipelineBase = 'PreCluster_Exp_ThresholdScan';    % name of CP pipeline that has been used to scan spot thresholds
S.ImportNameBase = 'PreSpots';                      % name given to the spots in the ScanSpotThresholds.m module

S.FolderCorrectionData = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\BiasCorrection');
S.strPlateVersion = 'Revision1';    % you can give a name to the scan that should be considered (note: in practice you might loop this code)

P{1}.name = '';     % optionally one might have traced a specific plate identifier in the Pipeline setting file and its outcome
P{1}.PlateDir = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleThresholdScan');

%%%%%%% SETTINGS FOR CLUSTER COMPUTING %%%%%%
strBrutus = 'bsub -W 8:00'; % standard submission command for LSF based clusters
strFunction = 'SpotThrDetection.brutusCorrectionOfPlateFromPipeline'; % function, which is called and actually obtains the correction matrix

%%%%%% PROCESSING %%%%%%%%%%
ensurePresenceOfDirectory(S.FolderCorrectionData);

currTime = datestr(clock,'yymmddHHMM');
S.currTime = currTime;
ClusterSettings.S = S;
ClusterSettings.P = P;

SettingsFileName = fullfile(...
    S.FolderCorrectionData ,['Submission_BiasCorrection_' ...
    currTime '.mat']);

save(SettingsFileName,'ClusterSettings');

numPlates = length(P);
for j=1:numPlates
    runDistributedJob(strBrutus, strFunction,SettingsFileName,j);
end


%%%%%%%%%%  one important detail !!!!! %%%%%%%%%%%%%%%
fprintf('Please note that, once the correction matrix is computed, you have to move the Correction Matrix to the BATCH folder of the plate, where it should be applied \n');