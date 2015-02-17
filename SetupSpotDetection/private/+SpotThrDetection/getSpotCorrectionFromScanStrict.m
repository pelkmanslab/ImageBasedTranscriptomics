function [CorrectionFunction ReferencePlaneRaw] = getSpotCorrectionFromScanStrict(ReferenceThreshold,ScannedThresholds,BatchPath,ImportName,ImageDimensions,SpotSmoothing,CorrectionSmoothing,scaleToInputThreshold)
% makes correction function for spot detection. It is derived from
% combining individual tested thresholds in such a way that the global spot
% count corresponds to the mean number of the initially selected Threshold
% calculation is done with single precision only to reduce memory load

% ReferenceThreshold = 0.055;              % Mean spots of this threshold will be used as reference for thresholds
% ScannedThresholds = 0.037:0.002:0.075;    % Thresholds that were tested [0.037:0.002:0.075]
% ImportName = 'ScanSpot';
% BatchPath = 'Y:\Data\Users\RNAFish\PLATES\120617_InSitu_CP102-1aa\BATCH';
% OutputName = 'CP102SpotCorr.mat';
% ImageDimensions = [2160 2560];

% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html

% default input
if nargin < 8
    scaleToInputThreshold = false;
else
    scaleToInputThreshold = logical(scaleToInputThreshold);
end

% Mildly filter (larger than spot size)
SmoothingSigma = 5;
H = fspecial('gaussian',[SpotSmoothing SpotSmoothing],SpotSmoothing./SmoothingSigma);

% Load Spot Number Per Pixel
numThresholds = length(ScannedThresholds);
smBias=zeros([ImageDimensions numThresholds],'single');

% [TS 120831] changed to single after checking that loss is negelectible to
% allow faster parallel computing of multiple thresholds on brutus without
% requesting more memory
% smBias=zeros([ImageDimensions numThresholds]);
% diffNumOcc = double(numOccurence) - single(numOccurence);
% diffSmBias = double(smBias) - single(smBias);

errorMessage = 'Please implement access to batch processed CP output, which follows the standards of your job manager. Currently only iBrain supported';

fprintf('Start to import location of spots.\n')
for k= 1:numThresholds
    ObjName  = [ImportName num2str(k)];
    strLocation = npc(fullfile(BatchPath,['Measurements_' ObjName '_Location.mat']));
if any(fileattrib(strLocation))   
    try
    numOccurence = VisualizePosition(strLocation,ObjName,1,false);
    catch 
        error(['Error during obtaining coordinates of spots. Potentially due to custom job manager.' errorMessage]);
    end
else
   fprintf(['Could not find ' strLocation '\n.']);
   error(errorMessage);
end
    smBias(:,:,k) = single(imfilter(numOccurence,H,'symmetric'));
    fprintf('%d out of %d Thresholds imported.\n',k,numThresholds)
end

fprintf('Start calculating the spot correction function.\n');

% Reference Value
[~, I] = min(abs(ScannedThresholds-ReferenceThreshold));
ReferencePlane = smBias(:,:,I);
ReferenceSpotNumber = mean(ReferencePlane(:));

% Reference Spots
ReferencePlane = I;


% Find Threshold closest to Reference
LeastDiff = abs(smBias-ReferenceSpotNumber);
[~, MinThr] = min(LeastDiff,[],3);

% Mildly blur image of thresholds
SmoothingSigma = 5;
H = fspecial('gaussian',[CorrectionSmoothing CorrectionSmoothing],CorrectionSmoothing*SmoothingSigma);
MinF =  imfilter(MinThr,H,'symmetric');

% interpolate
FlooredThresholdIndex = floor(MinF);
FlooredThresholdIndex(FlooredThresholdIndex==0) = 1;
FinalThresholds = ScannedThresholds(round(FlooredThresholdIndex)) + (MinF-FlooredThresholdIndex)*(ScannedThresholds(2)-ScannedThresholds(1));

% Create correction image. Note that filtered image will be divided by
% these values
if scaleToInputThreshold == false
    CorrectionFunction = double(FinalThresholds)./double(mean(FinalThresholds(:)));
else
    CorrectionFunction = double(FinalThresholds)./ ReferenceThreshold;
end

  ObjName  = [ImportName num2str(ReferencePlane)];
  strLocation = fullfile(BatchPath,['Measurements_' ObjName '_Location.mat']);
  ReferencePlaneRaw = double(VisualizePosition(strLocation,ObjName,1,false));

% 
% if bnShowFigure == true
%     
%     
%     RefQuant(1) = quantile(unCorrected(:),0.01);
%     RefQuant(2) = quantile(unCorrected(:),0.99);
%     
%     figure;
%     subplot(2,2,1);
%     imagesc(unCorrected, RefQuant);
%     title('Bias before correction');
%     colorbar;
%     subplot(2,2,2);
%     imagesc(CorrectionFunction, RefQuant);
%     title('Bias of correction function');
%     colorbar;
%     subplot(2,2,3);
%     hist(unCorrected(:),1000);
%     xlim(RefQuant);
%     title('Bias before correction');
%     subplot(2,2,4);
%     hist(CorrectionFunction(:),1000);
%     xlim(RefQuant);
%     title('Bias of correction function');
% end
end