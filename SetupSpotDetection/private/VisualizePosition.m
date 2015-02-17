function [numOccurence, numOccurenceScaled] = VisualizePosition(XorFile,YorMeasurmentName,downsampleFactor,bnNormalize,NrRandImagesPerOut)
% Visualizes the X and Y center of multiple objects


try
    if isnumeric(XorFile)
        YX = [YorMeasurmentName XorFile];
        clear YorMeasurmentName;
        clear XorFile;
    else
        load(XorFile);
        
        % sample random images
        if nargin<5
            doSampling = false;
        else
            if NrRandImagesPerOut <=length(handles.Measurements.(YorMeasurmentName).Location)
                doSampling = true;
            else
                doSampling = false;
            end
        end
        
        if doSampling == true;
            myRandomNumbers = rand(NrRandImagesPerOut,1);
            myRandomIndices = floor(myRandomNumbers*length(handles.Measurements.(YorMeasurmentName).Location));
            handles.Measurements.(YorMeasurmentName).Location = handles.Measurements.(YorMeasurmentName).Location(myRandomIndices);
            
        end
        
        
        
        %%%%%%
        
        
        
        numElements = cell2mat(cellfun(@(x) size(x,1), handles.Measurements.(YorMeasurmentName).Location, 'UniformOutput', 0));
        countPixels = sum(numElements,2);
        
        YX=NaN(countPixels,2);
        
        
        numLastIX=0;
        for k=1:size(handles.Measurements.(YorMeasurmentName).Location,2)
            
            numFirstIX=numLastIX+1;
            numLastIX=numLastIX+numElements(1,k);
            
            if ~isempty(handles.Measurements.(YorMeasurmentName).Location{1,k})
            [YX(numFirstIX:numLastIX,1)] = handles.Measurements.(YorMeasurmentName).Location{1,k}(:,2);
            [YX(numFirstIX:numLastIX,2)] = handles.Measurements.(YorMeasurmentName).Location{1,k}(:,1);
            end
            
            
        end
        clear handles;
    end
    
catch DataNotClear
    fprintf('Could not find input format /n')
end

YX = round(YX);

% remove data of centroids from image that is likely empty (CP output: [0 0]);
f  = (YX(:,1) == 0) & (YX(:,2) == 0);

YX = YX(~f,:);


% convert to linear index

dimYX = max(YX,[],1);

subYX=sub2ind2(dimYX,YX);
subYX = sort(subYX);


%%%%%%%



[uniqueSub, firstIX] = unique(subYX,'first');
lastIX = [firstIX(2:end);size(subYX,1)];
dOccurence = lastIX-firstIX+1;
clear subYX;

numOccurence = zeros(dimYX);

numOccurence(uniqueSub) = dOccurence;

if nargin<3
    downsampleFactor=10;
end
numOccurenceScaled = imresize(numOccurence, round(dimYX./downsampleFactor));

if nargin<4
    bnNormalize = false;
end


if bnNormalize == true
    numOccurenceScaled = numOccurenceScaled./mean(numOccurenceScaled(:));
end

minQuant = quantile(numOccurenceScaled(:),0.005);
maxQuant = quantile(numOccurenceScaled(:),0.995);
if nargout == 0
    imagesc(numOccurenceScaled, [minQuant maxQuant]);
    colorbar;
end

end