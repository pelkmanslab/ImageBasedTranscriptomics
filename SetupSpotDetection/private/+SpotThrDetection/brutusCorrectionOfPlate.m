function brutusCorrectionOfPlate(SettingsFile,PlateIndex)
load(SettingsFile);

P = ClusterSettings.P;
S = ClusterSettings.S;
k = PlateIndex;
ImportPathScan = ClusterSettings.ExportPath;
ExportPath = ImportPathScan;
P{k}.subsName

ReferenceThreshold = P{k}.numAimedThreshold;
TestedThresholds = P{k}.ThresholdRange;
ImportName = ['ScanSpot' P{k}.subsName];
OutputName = [P{k}.subsName 'SpotCorr'];
ImageDimensions = S.ImageDimensions;
SpotSmoothing = S.SpotSmoothing;
CorrectionSmoothing = S.CorrectionSmoothing;
bnShowFigure = false;

% Set Batch Path, either as Default (general impot + plate
% subfolder) or as BATCH of real plate

if isfield(P{k},'PlateBatchDir') == false
    BatchPath = fullfile(ImportPathScan,P{k}.SubfolderName);
else
    BatchPath = P{k}.PlateBatchDir;
end


OutObjectName = ['UnCorrected' P{k}.subsName];

OutputNameCorrection = fullfile(ExportPath, ['SpotCorrection' P{k}.subsName '.mat']);
OutputNameBiasAtAim = fullfile(ExportPath, ['SpotBiasUncorrected' P{k}.subsName '.mat']);
OutputNameSettings = fullfile(ExportPath, ['SpotDetectionSettings' P{k}.subsName '.mat']);

[CorrectionFunction unCorrected ReferencePlane] = SpotThrDetection.getSpotCorrectionFromScan(ReferenceThreshold,TestedThresholds,BatchPath,ImportName,OutputName,ImageDimensions,SpotSmoothing,CorrectionSmoothing,bnShowFigure);

SpotDetectionSettings.ReferencePlane = ReferencePlane;

SpotDetectionSettings.S = S;
SpotDetectionSettings.P = P;
SpotDetectionSettings.Id = PlateIndex;

save(OutputNameCorrection, 'CorrectionFunction');
save(OutputNameBiasAtAim,'unCorrected');
save(OutputNameSettings,'SpotDetectionSettings');





end