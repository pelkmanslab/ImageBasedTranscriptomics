function [CorrectionFunction unCorrected ReferencePlane] = getSpotCorrectionFromScan(ReferenceThreshold,TestedThresholds,BatchPath,ImportName,OutputName,ImageDimensions,SpotSmoothing,CorrectionSmoothing,bnShowFigure)
% makes correction function for spot detection. It is derived from
% combining individual tested thresholds in such a way that the global spot
% count corresponds to the mean number of the initially selected Threshold
% calculation is done with single precision only to reduce memory load

% ReferenceThreshold = 0.055;              % Mean spots of this threshold will be used as reference for thresholds
% TestedThresholds = 0.037:0.002:0.075;    % Thresholds that were tested [0.037:0.002:0.075]
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


% Mildly filter (larger than spot size)
if isempty(SpotSmoothing) || nargin<7;
    SmoothingSize = 20;
else
    SmoothingSize = SpotSmoothing;
end

SmoothingSigma = 5;
H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize./SmoothingSigma);

% Load Spot Number Per Pixel
numThresholds = length(TestedThresholds);
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
    strLocation = fullfile(BatchPath,['Measurements_' ObjName '_Location.mat']);
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
[~, I] = min(abs(TestedThresholds-ReferenceThreshold));
ReferencePlane = smBias(:,:,I);
ReferenceSpotNumber = mean(ReferencePlane(:));

% Reference Spots
ReferencePlane = I;


% Find Threshold closest to Reference
LeastDiff = abs(smBias-ReferenceSpotNumber);
[~, MinThr] = min(LeastDiff,[],3);

% Mildly blur image of thresholds
if isempty(CorrectionSmoothing) || nargin<8;
    SmoothingSize = 10;
else
    SmoothingSize = CorrectionSmoothing;
end
SmoothingSigma = 5;
H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);
MinF =  imfilter(MinThr,H,'symmetric');

% interpolate
FlooredThresholdIndex = floor(MinF);
FlooredThresholdIndex(FlooredThresholdIndex==0) = 1;
FinalThresholds = TestedThresholds(round(FlooredThresholdIndex)) + (MinF-FlooredThresholdIndex)*(TestedThresholds(2)-TestedThresholds(1));

% Create correction image. Note that filtered image will be divided by
% these values
CorrectionFunction = double(FinalThresholds)./double(mean(FinalThresholds(:)));

if strcmp(OutputName((end-3):end),'.mat') == false
    OutputName = [OutputName '.mat'];
end

try
    save(fullfile(BatchPath,OutputName),'CorrectionFunction');
    fprintf(['Saved' OutputName '.\n']);
catch notSave
    fprintf(['Could not save' OutputName ' Note that output has to be .mat and target destination writeable.\n']);
end


% Output Figure
if nargin<9;
    bnShowFigure = true;
end

unCorrected = double(ReferencePlane)./double(ReferenceSpotNumber);

ReferencePlane = double(ReferencePlane);

if bnShowFigure == true
    
    
    RefQuant(1) = quantile(unCorrected(:),0.01);
    RefQuant(2) = quantile(unCorrected(:),0.99);
    
    figure;
    subplot(2,2,1);
    imagesc(unCorrected, RefQuant);
    title('Bias before correction');
    colorbar;
    subplot(2,2,2);
    imagesc(CorrectionFunction, RefQuant);
    title('Bias of correction function');
    colorbar;
    subplot(2,2,3);
    hist(unCorrected(:),1000);
    xlim(RefQuant);
    title('Bias before correction');
    subplot(2,2,4);
    hist(CorrectionFunction(:),1000);
    xlim(RefQuant);
    title('Bias of correction function');
end
end