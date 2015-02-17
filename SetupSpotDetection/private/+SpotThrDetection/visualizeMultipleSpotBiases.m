function visualizeMultipleSpotBiases(cellBatchPath,cellObjName,cellDescription,ImageDimensions)
% visualizes Spot Bias of multiple Spot Detections Note that quantiles of
% first File will be used for scaling


% cellBatchPath{1} = 'W:\Data\Users\RNAFish\PLATES\120617_InSitu_CP102-1aa\Previous Pipelines\PreAnalysis 02 RNAFishPipeline_120707\BATCH';
% cellObjName{1} = 'SpotsM';
% cellDescription{1} = 'CP102-1aa before Correction with flexible Thresholds';
% 
% cellBatchPath{2} = 'W:\Data\Users\RNAFish\PLATES\120617_InSitu_CP102-1aa\BATCH';
% cellObjName{2} = 'DSpots';
% cellDescription{2} = 'CP102-1aa after Correction with flexible Thresholds';
%
% ImageDimensions = [2160 2560];      % dimensions of imagec



% Smoth, note that quite large to see overall bias
% SmoothingSize = 100; % changed by [TS 130219] to make larger trends
%visible
SmoothingSize = 150;


SmoothingSigma = 5;
H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize./SmoothingSigma);

% Load Spot Number Per Pixel
numSpotfiles = length(cellBatchPath);
smBias=zeros([ImageDimensions numSpotfiles],'single');
fprintf('Start to import location of spots.\n')
for k= 1:numSpotfiles
    strLocation = fullfile(cellBatchPath{k},['Measurements_' cellObjName{k} '_Location.mat']);
     numOccurence = VisualizePosition(strLocation,cellObjName{k},1,false);
%     numOccurence = CountObjectsPerPixel(strLocation,'cv7k');
    filtImage = imfilter(numOccurence,H,'symmetric');
    smBias(:,:,k) = single(filtImage./mean(filtImage(:)));
    fprintf('%d out of %d Spot Detections imported.\n',k,numSpotfiles)       
end

% Output Figure
firstImage = smBias(:,:,1);

RefQuant(1) = quantile(firstImage(:),0.01);
RefQuant(2) = quantile(firstImage(:),0.99);

figure;
for k=1:numSpotfiles
    currImage = smBias(:,:,k);
    
    subplot(2,numSpotfiles,k);
    imagesc(currImage, RefQuant);
    title([cellObjName{k} ' : ' cellDescription{k}]);
    colorbar;

    subplot(2,numSpotfiles,(numSpotfiles+k));
    hist(currImage(:),1000);
    xlim(RefQuant);
    title([cellObjName{k} ' : ' cellDescription{k}]);
    xlabel('Spot Count / mean Spot Count')
    ylabel('Pixels');
end

end