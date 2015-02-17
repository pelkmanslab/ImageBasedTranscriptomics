function brutusCorrectionOfPlateFromPipeline(SettingsFile,PlateIndex)
load(SettingsFile);


P = ClusterSettings.P;
S = ClusterSettings.S;
k = PlateIndex;
if ischar(k)
   k = str2num(k); 
end

CurrName = P{k}.name

% check input
if ~isfield(S,'scaleToInputThreshold')
    scaleToInputThreshold = false;
else
    scaleToInputThreshold = S.scaleToInputThreshold;
end


% load pipeline
strPlateDir = P{k}.PlateDir;
strPipeline = nnpc(fullfile(strPlateDir,[S.PipelineBase CurrName '.mat']));
if any(strPipeline)
    Precluster = load(strPipeline);
else
   strPipeline
   error('could not fine pipeline');
end


% get information from pipeline
IXScanModule = find(cell2mat(cellfun(@(x) strcmp(x,'ScanSpotThresholds'), Precluster.Settings.ModuleNames, 'UniformOutput', false)));

FilterSize =            str2num(Precluster.Settings.VariableValues{IXScanModule,4})
MinMaxQuant =           str2num(Precluster.Settings.VariableValues{IXScanModule,5})
RescaleLimits =         str2num(Precluster.Settings.VariableValues{IXScanModule,6})
ScannedThresholds =     str2num(Precluster.Settings.VariableValues{IXScanModule,7})
PixelIntensityLimit =   str2num(Precluster.Settings.VariableValues{IXScanModule,9})
ReferenceThreshold =    str2num(Precluster.Settings.VariableValues{IXScanModule,12})


% settings for spot bias correction
ImageDimensions = S.ImageDimensions;
SpotSmoothing = S.SpotSmoothing;
CorrectionSmoothing = S.CorrectionSmoothing;


BatchPath = fullfile(strPlateDir,'BATCH');
ImpName = [S.ImportNameBase CurrName];



[CorrectionFunction ReferencePlaneRaw] = SpotThrDetection.getSpotCorrectionFromScanStrict(...
    ReferenceThreshold,ScannedThresholds,BatchPath,ImpName,ImageDimensions,SpotSmoothing,CorrectionSmoothing,scaleToInputThreshold);


RefSettings.FilterSize = FilterSize;
RefSettings.MinMaxQuant = MinMaxQuant;
RefSettings.RescaleLimits = RescaleLimits;
RefSettings.ScannedThresholds = ScannedThresholds;
RefSettings.PixelIntensityLimit = PixelIntensityLimit;
RefSettings.ReferenceThreshold = ReferenceThreshold;
RefSettings.ClusterSettings = ClusterSettings;
RefSettings.CurrJob = k;
RefSettings.CurrName = CurrName;




ExportPath = S.FolderCorrectionData;

OutputnameCorrection = fullfile(ExportPath, ['SpotBiasCorrection_' CurrName '.mat']);
OutputnameReferencePlane = fullfile(ExportPath, ['ReferencePlane_' CurrName '.mat']);
OutputnameParameters = fullfile(ExportPath, ['ParametersOfScan_' CurrName '.mat']);



save(OutputnameCorrection, 'CorrectionFunction');
save(OutputnameReferencePlane,'ReferencePlaneRaw');
save(OutputnameParameters,'RefSettings');





end