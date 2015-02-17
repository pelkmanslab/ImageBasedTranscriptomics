function [ObjCount] = runObjByFilterOnMultipleImages(ImageNames,PlatePath,Filter,ObjThr,limQuant,RescaleThr,ObjIntensityThr,closeHoles,ObjSizeThr,IllumCorrRef, fractionImages,DetectionBias)
% will run Object on the images provided in the cell IMAGENAMES, which
% has names of images of one Image Set grouped togehter. Output will be
% CELL  with matrices. Where rows denote different sites/images and columns
% denote the Object Count for each value provided by ObjThr


% if no Detection Bias is set, set to empty
if nargin<12
    DetectionBias = [];
end

numImageGroups = size(ImageNames,2);
PresentImageSets = cell2mat(cellfun(@(x) size(x,1), ImageNames, 'UniformOutput', false));
NumImagesToAnalye = ceil(PresentImageSets.*fractionImages);

% Determine if 2D or 3D images
if size(ImageNames{1},2) > 1
    imageMode = 'CP3D';
else
    imageMode = 'CPClassical';
end

% Check if illumination correction should be applied
if ~isempty(IllumCorrRef.IllumFilt_STD)
    bnApplyIlluminationCorrection = true;
else
    bnApplyIlluminationCorrection = false;
end

% Determine number of objects
ObjCount = cell(numImageGroups,1);
for l=1:numImageGroups
    ObjCount{l} = nan(NumImagesToAnalye(l),length(ObjThr));
    randIX = randperm(size(ImageNames{l},1));
    for k = 1:NumImagesToAnalye(l)
        
        switch imageMode
            case 'CPClassical'
                cFile = ImageNames{l}(randIX(k),:);
                strImage = fullfile(PlatePath,'/TIFF/',cFile{:} );
            case 'CP3D'
                strImage = cellfun(@(x) fullfile(PlatePath,'/TIFF',x), ImageNames{l}(randIX(k),:),'UniformOutput', false);
        end
         
        % load Image
        try
            switch imageMode
                case 'CP3D'
                    strImage = cellfun(@(x) nnpc(x), strImage, 'UniformOutput', false); % conversion of path useful for local debugging
                    OrigImage = imreadCP3D(strImage);
                case 'CPClassical'
                    OrigImage = imread(strImage);
            end
        catch notLoaded
            error(['Could not load Image with file path/name ' strImage '.'])
        end
        
        
        % apply Illumination correction, if requested
        if bnApplyIlluminationCorrection == true
            
            switch size(OrigImage,3)
                case 1    % in case of 2D image use Nico's function for illum correction
                    OrigImage = IllumCorrect(OrigImage,IllumCorrRef.IllumFilt_Mean,IllumCorrRef.IllumFilt_STD,1);
                otherwise % otherwise use implementation supporting 3D
                    OrigImage = applyNBBSIllumCorrCP3D(OrigImage,IllumCorrRef.IllumFilt_Mean,IllumCorrRef.IllumFilt_STD);
            end
        end
        
        % Get object counts at different thresholds
        [ObjCount{l}(k,:)] = ObjByFilter(OrigImage,Filter,ObjThr,limQuant,RescaleThr,ObjIntensityThr,closeHoles,ObjSizeThr,DetectionBias);
        
        
        FractionAnalyzsed = k./NumImagesToAnalye(l).*100;
        if mod(k,floor(NumImagesToAnalye(l)./10))==0
            fprintf('%0.0f%% of image set %d of %d \n',FractionAnalyzsed,l,numImageGroups);
        end
    end
end

end