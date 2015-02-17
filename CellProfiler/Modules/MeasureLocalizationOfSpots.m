function handles = MeasureLocalizationOfSpots(handles)

% Help for the Measure Localization Of Spots Features module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures the spot localization features in the cytoplasm of
% cells as destribed by Battich et al., 2013.
% *************************************************************************
%
% Requirement: Children objects (such as Spots) have to be "related" to
% parents (such as Cells) by a Relate module, prior to this module.
%
% Given the children objects (e.g. cytoplasmic spots), the parent objects
% (e.g. cells) and a third reference object (e.g. nuclei), this module
% computes the localization features of each children with respect with its
% parent, the third object provided and every other children object of the
% parent.
%
% Measurement:                                              Feature Number:
% Closest distance to parent outline                       |       1
% Neighbour status of the closest membrane
% (1 is it close to a neighbour membrane and
% 2 it is close to a non neighbour membrane, if the
% border of the image is closest this is set to nan)       |       2
% Distance of children to parent Centroid                  |       3
% Distance of children to third reference Centroid         |       4
% Distance of children to membrane along
% the children-third reference axis                        |       5
% Mean distance of child to all other children in parent   |       6
% Std of dist. of child to all other children in parent    |       7
% Variance of dist. of child to all other children in par. |       8
% Distances from child to include x number of percentages
% of other children in the parent                          |     9->8+x
% Fraction of all other children at y distances from child | 8+x+1->8+x+y
%
%
% Note: all distances are calculated in pixels. 
%
% Two sets of outputs are generated:
% x raw Measurments
% x Measurmetents, where data of each spot is normalized by Z-scoring
% against 100 iteratations of randomly chosen pixels, which are not in the
% tertiary object
%
% ADJUSTMENTS for LOW THROUGHPUT platforms:
% This module has been designed for large scale image based
% transcriptomics. Lower throughput analysis might require two changes
% depending on microscopic / computational infrastructure. You can find
% these position in the code by searching for 
%    1) Exclusion of false cells
%    and
%    2) Number of randomizations
%
% ad 1): If a cell is larger than half of the image it is excluded. On
% large field of view cameras, there is not a single real cell which is as
% large (as screening setups will include approx. 50-200 cells per image).
% The only "cell" objects of that size, which we observed, were very rare, 
% falsely identified objects / background. Running this module will
% increase memory requiremt by around 30GB for such objects (causing 
% unnecessary resubmission to more powerful computational nodes / inability
% to continue calculation).
% With high magnification lenses and small field of view cameras,
% biological cells might be larger than half of the image. Thus the limit
% has to be changed or removed.
%
% ad 2): while processing on a single computer will last less than 3min for
% the majority of all high abundance genes, extremely high abundance RNAs,
% such as actin or ribosomal RNA, will take around 20min. The overall
% processing speed of this module can be increased massively by reducing
% the number of random iterizations, which are used for spot normalization.
% Be warned however, that this will result in a less accurate estimate of the
% measurements and should usually be avoided. If you are only intrested in
% raw measurments, you might decrease the number of random iterations to 1,
% which will speed up the time of the module by a factor of roughly 100
% (excluding 1 min of preprocessing)
%
%
%
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%

%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What are the parent objects? e.g Cells.
%infotypeVAR01 = objectgroup
%defaultVAR01 = Cells
%choiceVAR01 =
ParentObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What are the children objects for which you what to compute distances?
%infotypeVAR02 = objectgroup
%defaultVAR02 = Spots
%choiceVAR02 =
ChildrenObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What are the objects for which you what to compute distance to the children objects? (note: this object must be a direct child or a direct parent of the 'Parent' inputed above, e.g. Nuclei)
%infotypeVAR03 = objectgroup
%defaultVAR03 = Nuclei
%choiceVAR03 =
ThirdObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What is the distance between two neighbouring cells? (pixel) recommended: 10
%defaultVAR04 = 10
DilationFactor = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Distance(s) required to cover this / these fraction(s) of sibling spots (values between 0 and 1).
%defaultVAR05 = [0.10 0.25 0.50 0.75]
PercentageVect = (handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Fraction of siblings at this / these distances (pixels).
%defaultVAR06 = [40 90]
RadiousVect = (handles.Settings.VariableValues{CurrentModuleNum,6});



%%%%%%%%%%%%%%%%%%%%
%%% Inputs check %%%
%%%%%%%%%%%%%%%%%%%%

%%% check the parent object
if isempty(ParentObjectName)
    error('%s: Please enter a Parent object.',mfilename);
end

%%% check the child object
if isempty(ChildrenObjectName)
    error('%s: Please enter a Child object.',mfilename);
end

%%% check for any third refernce object, if abscent give a warning
ThirdStatus = true;
if isempty(ThirdObjectName)
    warning('%s: No third object entered.',mfilename);
    ThirdStatus = false;
end

%%% check if children have been related to the parent objects
if isfield(handles.Measurements.(ChildrenObjectName),'Parent')
    IXParentObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.(ChildrenObjectName).ParentFeatures, 'uniformoutput',false)),1,'first');
    if isempty(IXParentObject)
        error('%s: Please relate the clidren objects to the parent objects.',mfilename)
    end
else
    error('%s: Please relate the clidren objects to the parent objects.',mfilename)
end


%%% evaluate radious and percentage vectors and check that no non-allowed
%%% charecter has been introducaed
[isSafe PercentageVect]= inputVectorsForEvalCP3D(PercentageVect,false);
if isSafe==false
    error(['Please enter vector in correct format'])
end
matPerVector = eval(PercentageVect);

[isSafe RadiousVect]= inputVectorsForEvalCP3D(RadiousVect,false);
if isSafe ==false
    error(['Please enter vector in correct format'])
end
matRadiousVector = eval(RadiousVect);


%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization %%%
%%%%%%%%%%%%%%%%%%%%%%


% create names of features.
if ThirdStatus
    
    strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
        strcat('Label_Outline_',ParentObjectName),...
        strcat('Dis_Centroid_',ParentObjectName),...
        strcat('Dis_Centroid_',ThirdObjectName),...
        strcat('Dis_Proj_',ThirdObjectName,'_',ParentObjectName),'Dis_Mean','Dis_Std','Dis_Var'};
    for i = 1:length(matPerVector)
        strFeautesGlobal{end+1} = strcat('Dis_DistForFraction_',sprintf('%.2f',matPerVector(i)));
    end
    for i = 1:length(matRadiousVector)
        strFeautesGlobal{end+1} = strcat('Dis_FractionAtDist_',sprintf('%.4d',matRadiousVector(i)));
    end
    
    strFeautesGlobalRaw = cellfun(@(x) [x 'Raw'], strFeautesGlobal,'UniformOutput',false);
    strFeautesGlobalZScored = cellfun(@(x) [x 'ZScored'], strFeautesGlobal,'UniformOutput',false);
    
else
    
    strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
        strcat('Label_Outline_',ParentObjectName),...
        strcat('Dis_Centroid_',ParentObjectName),'Dis_Mean','Dis_Std','Dis_Var'};
    for i = 1:length(matPerVector)
        strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.2f',matPerVector(i))); %#ok<*AGROW>
    end
    for i = 1:length(matRadiousVector)
        strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.4d',matRadiousVector(i)));
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing of the Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% retrieve segmentation images
SegmentedParentObjectImage = CPretrieveimage(handles,['Segmented', ParentObjectName],ModuleName,'MustBeGray','DontCheckScale');
SegmentedChildrenObjectImage = CPretrieveimage(handles,['Segmented', ChildrenObjectName],ModuleName,'MustBeGray','DontCheckScale');
SegmentedThirdObjectImage = CPretrieveimage(handles,['Segmented', ThirdObjectName],ModuleName,'MustBeGray','DontCheckScale');


%%% check if the image is empty
IXObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.Image.ObjectCountFeatures, 'uniformoutput',false)),1,'first');
ObjCounts = handles.Measurements.Image.ObjectCount{1,handles.Current.SetBeingAnalyzed}(1,IXObject);
ImageEmpty = false;

if ObjCounts < 1
    ImageEmpty = true;
    listUniq = [];
else
    listUniq = [1:ObjCounts];
end


%%% initialize parameters
%%% retrieve the location of the children
matChildrenLocation = handles.Measurements.(ChildrenObjectName).Location{handles.Current.SetBeingAnalyzed};

SizeOfChildren = size(matChildrenLocation,1);
numTotalFeatures = length(strFeautesGlobal);
matLocMeasurements = NaN(SizeOfChildren,numTotalFeatures);
matLocMeasurementsMean = NaN(SizeOfChildren,numTotalFeatures);
matLocMeasurementsStd = NaN(SizeOfChildren,numTotalFeatures);
matLocMeasurementsZScored = NaN(SizeOfChildren,numTotalFeatures);


clear SizeOfChildren



%%% if there are any parent and children objects proceed to calculate
%%% distance features

if ImageEmpty == false
    
    % Obtain the coordinates of pixels within each parent
    ParentPixels = regionprops(SegmentedParentObjectImage,'PixelList');
    SegmentedNonThird = SegmentedParentObjectImage.*(SegmentedThirdObjectImage<=0);
    NonThirdPixels = regionprops(SegmentedNonThird,'PixelList');
    
    
    fprintf('%s: calculating the final border image\n',mfilename)
    
    %%% check background. It was reported that sometimes the background is not 0
    listUniqInImage = unique(SegmentedParentObjectImage);
    backgroundID = listUniqInImage(~ismember(listUniqInImage',listUniq)');
    
    
    %%% do default expansion of cells for neighbour calculation
    [SegmentedParentObjectImageExp,CurrentObjLabels] = bwdist(SegmentedParentObjectImage);
    SegmentedParentObjectImageExp = (SegmentedParentObjectImageExp < DilationFactor).*SegmentedParentObjectImage(CurrentObjLabels);
    
    isBackground=false(size(SegmentedParentObjectImageExp));
    for k=1:length(backgroundID)
        isBackground(SegmentedParentObjectImageExp==backgroundID(k))=true;
    end
    SegmentedParentObjectImageExp(isBackground)=0;
    
    
    %%% expand background by two pixels
    isBackgroundExp = bwmorph(isBackground,'dilate',2);
    
    %%% get the edges of the cells, in such a way that later shrinking
    %%% would still yield a continuous edge
    filter1 = fspecial('sobel'); %
    filter2 = filter1';
    %%% Applies each of the sobel filters to the original image.
    I1 = imfilter(SegmentedParentObjectImageExp, filter1);
    I2 = imfilter(SegmentedParentObjectImageExp, filter2);
    %%% Adds the two images. Use it for edge detection
    matEdgeImagesIX = (abs(I1) + abs(I2))>0; clear I1; clear I2;
    
    %%% get edge neighbour and edge non-neighbour info
    matEdge_NonEdgeIX = and(isBackgroundExp,matEdgeImagesIX);
    matEdge_NonEdge = double(matEdgeImagesIX);
    matEdge_NonEdge(matEdge_NonEdgeIX) = 2;
    matEdge = bwmorph(matEdge_NonEdge==1,'dilate',1);
    matNonEdge = bwmorph(matEdge_NonEdge==2,'dilate',1);
    matEdge_NonEdge(matEdge) = 1;
    matEdge_NonEdge(matNonEdge) = 2;
    
    %%% set all borders of an image to the value of 3 ; if using double
    %%% sobel: instead of only one pixel: 2
    matEdge_NonEdge(:,1:2) = 3;
    matEdge_NonEdge(:,(end-1:end)) = 3;
    matEdge_NonEdge(1:2,:) = 3;
    matEdge_NonEdge((end-1:end),:) = 3;
    
    %     figure;imagesc(matEdge_NonEdge) Debug
    
    %%% calculate allowed coordinates per object (to reduce computational cost)
    %-> recalculating objects is necesary because of a background bug!
    LisOfObjects = unique(SegmentedParentObjectImageExp(:));
    if LisOfObjects(1) == 0 % remove background, if present
        LisOfObjects = LisOfObjects(2:end);
    end
        
    
    %%% Obtain pixels at inner periphery of cells. note that this will
    %%% reduce the computing time of spot features by more than 50% since much
    %%% less membrane pixels have to be considered
    distanceToObjectMax = 3;
    FinalLabelMatrixImage = SegmentedParentObjectImage;
    props = regionprops(FinalLabelMatrixImage,'BoundingBox');
    BoxPerObj = cat(1,props.BoundingBox);
    
    % get outer coordinates of bounding box of each object (note that
    % bounding boxes will speed up morphological image operations)
    N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                                 N(f) = 1;
    S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(FinalLabelMatrixImage,1);    S(f) = size(FinalLabelMatrixImage,1);
    W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                                W(f) = 1;
    E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(FinalLabelMatrixImage,2);    E(f) = size(FinalLabelMatrixImage,2);
    
    % create empty output
    FinalLabelMatrixImage2  = zeros(size(FinalLabelMatrixImage));
    if ~isempty(LisOfObjects)  % if objects present
        clear k;
        for k=LisOfObjects'  % loop through individual objects to safe computation
            
            % create small image to speed up morphological operations
            miniImage = FinalLabelMatrixImage(N(k):S(k),W(k):E(k));
            bwminiImage = miniImage==k;
            
            % get pixels within object which are next to background. note
            % that the alternative strategy to shrink the object could lead
            % to small parts of the membrane becoming disconnected, if
            % there is a finger like extension of cells. When following the
            % reverse strategy to extend the background there can be a few
            % extra pixels, however since they are rare they will be
            % neglectible for the speed of the later calculation and it is
            % favourable to have closed cell membranes
            padbw = padarray(bwminiImage,[1 1]); % pad with 0s to ensure that background on all sides
            ipadbw = ~padbw; % get background
            tpadbw = bwmorph(ipadbw,'dilate',1); % extend background
            dpadbw = tpadbw & ~ipadbw; % get inner pixels of objects
            dpadbwR = dpadbw(2:(end-1),2:(end-1)); % reverse padding to get correct coordinates
            
            dpadbwR = bwmorph(dpadbwR,'thin'); % thin so that across diagonals around 1/3 of pixels is lost. will increase speed spot parameter calculation.
            
            % now map back the linear indices
            [r c] = find(dpadbwR);
            
            % get indices for final image (note that mini image might have
            % permitted regions of other cells and thus boxes can not be 
            % directly overlaid).
            r = r-1+N(k);
            c = c-1+W(k);
            w = sub2ind(size(FinalLabelMatrixImage2),r,c);
            
            % Update Working copy of Final Segmentation image based on
            % linear indices. (Linear indices prevent depenceny of other
            % objects in same box. Thus overlapping boxes
            % can be processed sequentially)
            FinalLabelMatrixImage2(w) = k;
            
        end
    end
    
    cellMembraneLocationId = regionprops(FinalLabelMatrixImage2,'PixelList');
    
    clear FinalLabelMatrixImage2; clear FinalLabelMatrixImage;
    
    cellMembraneLocation = cell(LisOfObjects(end),1);
    for j=LisOfObjects'
        CurrCand = cellMembraneLocationId(j).PixelList;
        if ~isempty(CurrCand)
            cellMembraneLocation(j) = {[CurrCand(:,2) CurrCand(:,1)]}; % preserve format of previous versions of this module and use [Y X] for membrane coordinates};
        end
    end
    
    % Update Edge image so that it is specified for each pixel (note that
    % membrane pixels are now within cells); While calucating label for
    % each pixel, which depends upon the status of the closest edge, is
    % slow, it is only done once per image. It allows to use inner pixels
    % of the cells as membrane pixels - even if they might be away from
    % edge in extended segementation image by a few pixels. This way the
    % size of the cells does not have to become changed and distance does
    % not have to caluclated back.
    [~,ELabel] = bwdist(matEdge_NonEdge);
    matEdge_NonEdge = matEdge_NonEdge(ELabel);
    clear ELabel;
    
    % label pixels closest to image border with nan instad of 3. Thus only
    % two values will be used as a label. Thus nan robust meausrements will
    % be directly interpretable
    matEdge_NonEdge(matEdge_NonEdge==3) = nan;
    
    %get the ID of the parents of each child (spot)
    matParentObject = handles.Measurements.(ChildrenObjectName).Parent{handles.Current.SetBeingAnalyzed}(:,IXParentObject);
    
    %%% get parent  and third object centroid coordinates
    matParentLocation = handles.Measurements.(ParentObjectName).Location{handles.Current.SetBeingAnalyzed};
    matThirdLocation = handles.Measurements.(ThirdObjectName).Location{handles.Current.SetBeingAnalyzed};
    
    %%% Calculate distances
    fprintf('s: Calculating first set of distances\n',mfilename)
    
    %Include safety check for cells, that must not be analysed
    
    %%% A) The background
    forbiddenObjectIDs = 0;
    
    %%% B) DISCARD CELLS > 50% of image
    %    1) Exclusion of false cells
    % Cells, which occupy more than half of the image are discarded. In our
    % case these represent rare cases, where cells are not correctly
    % identified. They will ususally increase the memory requirement by
    % >30GB, which causes some computational jobs to crash. If these cells
    % should not be discarded, consider commenting the following code or
    % changing the limit
    nonZeroIds = sort(SegmentedParentObjectImageExp(SegmentedParentObjectImageExp~=0)); % get non background pixels and sort them
    if any(nonZeroIds)
        
        tmpAreaOfCells = regionprops(SegmentedParentObjectImageExp,'Area');
        totpixel = length(SegmentedParentObjectImage(:));
        for jj =1: length(tmpAreaOfCells)
            if(tmpAreaOfCells(jj).Area/totpixel) > 0.5  % 1) Exclusion of false cells: replace 0.5 by 1 to prevent any discarding
                forbiddenObjectIDs(end+1) = jj; %#ok<AGROW>
                fprintf('At least one parent object is larger than half of the image. This Objcect will be excluded from analysis. \n')

            end
        end
        clear tmpAreaOfCells;
        
    end
    clear nonZeroIds; clear SegmentedParentObjectImageExp;
    % end of B) DISCARD CELLS > 50% image
    
    tic % on actin images, with many spots, calculation should be around 20 min
    %%% looping through cells
    for i = unique(matParentObject)'
        
        %%% only do calculations if the children do not fall in the
        %%% forbidden and if there is no bug in segmentation
        if ~any(ismember(i,forbiddenObjectIDs)) && ~isempty(cellMembraneLocation{i})
            
            %get all the children of the parent
            tempChildIx = matParentObject == i;
            tempTotalChildren = sum(tempChildIx);
            tempParentLocation = matParentLocation(i,:);
            tempChildrenLocation = matChildrenLocation(tempChildIx,:);
            
            if tempTotalChildren > 0
                %%% Real spots %%%
                CurrCellMembraneLocalization = cellMembraneLocation{i};
                tempThirdLoc = matThirdLocation(i,:);
                
                % Obtain Measurments
                matLocMeasurementsSingleCellSingleIteration = CPgetSpotLocalizations(tempParentLocation,tempChildrenLocation,tempThirdLoc,...
                    CurrCellMembraneLocalization,tempTotalChildren,matEdge_NonEdge,matPerVector,matRadiousVector,i);
                matLocMeasurements(tempChildIx,:) = matLocMeasurementsSingleCellSingleIteration;
                
                %%% Standardization %%%
                % Use radomly placed pixels as spots
                numIterations = 100; %2) Number of randomizations : set the number of random iteratizations. Might adjust according to your computational infrastructure.
                CurrParentsPixels =  ParentPixels(i).PixelList;
                if size(CurrParentsPixels,1) < tempTotalChildren % if the region, which is not tertiary object is smaller than the amount of children, use pixels of complete cell as reference
                    CurrParentsPixels = NonThirdPixels(i).PixelList;
                end
                
                
                numParentPixels = size(CurrParentsPixels,1);
                matLocMeasurementsSingleCellRandomSpots = NaN(tempTotalChildren,numTotalFeatures,numIterations);
                
                for jj=1:numIterations
                    % Get random pixel. The amount will correspond to the
                    % amount of real spots. As with real spots they are not
                    % overlapping. Note about algorithm: permutate small
                    % amount of pixels first and hope that they are non
                    % overlapping
                    randIX = ceil(max(numParentPixels)*rand(max(ceil(tempTotalChildren*1.2)),1));
                    [u IX] = unique(randIX);
                    if length(u)>= tempTotalChildren
                        IX = IX(randperm(length(IX)));
                        randIX = randIX(IX(1:tempTotalChildren));
                    else % If not, permute all pixels. Note that this would be slow, if executed for all cells
                        randIX = randperm(numParentPixels); 
                        randIX  = randIX(1:tempTotalChildren);
                    end
                    tempRandomChildLocation = CurrParentsPixels(randIX,:);
                    
                    % Obtain Measurements. note that they are only kept in
                    % memory for each individual cell
                    matLocMeasurementsSingleCellRandomSpots(:,:,jj) = CPgetSpotLocalizations(tempParentLocation,tempRandomChildLocation,tempThirdLoc,...
                        CurrCellMembraneLocalization,tempTotalChildren,matEdge_NonEdge,matPerVector,matRadiousVector);
                end
                
            end
            
            % Z score for each measurment of each individual really
            % measured spot. In principle there would be two ways for the
            % standardization a) against one spot of each repetition. b)
            % against all spots of the cell of each repetition. ; While in
            % b) the measurments might be more accurate, it would mean that
            % secondary features, such as the variability of the spot
            % features would indirectly mirror the amount of spots, if the
            % amount of spots is different. Thus a) is used
            matLocMeasurementsMean(tempChildIx,:) = nanmean(matLocMeasurementsSingleCellRandomSpots,3);
            matLocMeasurementsStd(tempChildIx,:) = nanstd(matLocMeasurementsSingleCellRandomSpots,[],3);
            
            
        end
    end
    toc
    
    matLocMeasurementsZScored = (matLocMeasurements-matLocMeasurementsMean)./matLocMeasurementsStd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% use format of last version of CellProfiler 1
%
%       """""""""""  WARNING: the following lines are outdated and have to
%       be adjusted for to use columsn of matLocMeasurements and
%       matLocMeasurementsZScored since individual outputs are now saved
%      in one large matrix by the subfunction CPgetSpotLocalizations """""
%      also they do not contain normalized data
%
% if ThirdStatus : note
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{1}]){SetBeingAnalyzed} = matDistancesToOutline;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{2}]){SetBeingAnalyzed} = matOutlineLabel;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{3}]){SetBeingAnalyzed} = matDisToParentCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{4}]){SetBeingAnalyzed} = matDisToThirdCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{5}]){SetBeingAnalyzed} = matDisToOutlineThirdLine;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6}]){SetBeingAnalyzed} = matMeanDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{7}]){SetBeingAnalyzed} = matStdDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8}]){SetBeingAnalyzed} = matVarDis;
%     for i = 1:length(matPerVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8+i}]){SetBeingAnalyzed} = matDisPerX1(:,i);
%     end
%     for i = 1:length(matRadiousVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8+length(matPerVector)+i}]){SetBeingAnalyzed} = matDisRadY1(:,i);
%     end
% else
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{1}]){SetBeingAnalyzed} = matDistancesToOutline;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{2}]){SetBeingAnalyzed} = matOutlineLabel;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{3}]){SetBeingAnalyzed} = matDisToParentCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{4}]){SetBeingAnalyzed} = matMeanDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{5}]){SetBeingAnalyzed} = matStdDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6}]){SetBeingAnalyzed} = matVarDis;
%     for i = 1:length(matPerVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6+i}]){SetBeingAnalyzed} = matDisPerX1(:,i);
%     end
%     for i = 1:length(matRadiousVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6+length(matPerVector)+i}]){SetBeingAnalyzed} = matDisRadY1(:,i);
%     end
% end



%%% if older versions of CellProfiler are being used, you can save the
%%% measurments in this way

%%% save the features names
if handles.Current.SetBeingAnalyzed==1
    
    handles.Measurements.(ChildrenObjectName).ChildLocalizationRawFeatures = strFeautesGlobalRaw;
    handles.Measurements.(ChildrenObjectName).ChildLocalizationRaw = cell(1,handles.Current.NumberOfImageSets);
    
    handles.Measurements.(ChildrenObjectName).ChildLocalizationZScoredFeatures = strFeautesGlobalZScored;
    handles.Measurements.(ChildrenObjectName).ChildLocalizationZScored = cell(1,handles.Current.NumberOfImageSets);
    
end

%%% save measurements
%   commented out: code from previous version
%     handles.Measurements.(ChildrenObjectName).ChildrenLocalization{handles.Current.SetBeingAnalyzed} = ...
%         [matDistancesToOutline ...
%         matOutlineLabel ...
%         matDisToParentCentroid ...
%         matDisToThirdCentroid ...
%         matDisToOutlineThirdLine ...
%         matMeanDis ...
%         matStdDis ...
%         matVarDis ...
%         matDisPerX1 ...
%         matDisRadY1];
% handles.Measurements.(ChildrenObjectName).ChildrenLocalization{handles.Cu
% rrent.SetBeingAnalyzed} = matLocMeasurementsZScored; % commetned out
% saving of previous version. note that name of features is changed so that
% modules can be run in parallel
handles.Measurements.(ChildrenObjectName).ChildLocalizationRaw{handles.Current.SetBeingAnalyzed} = matLocMeasurements;
handles.Measurements.(ChildrenObjectName).ChildLocalizationZScored{handles.Current.SetBeingAnalyzed} = matLocMeasurementsZScored;


try
    %Visuallize outputs
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        CPresizefigure(SegmentedParentObjectImage,'TwoByTwo',ThisModuleFigureNumber);
        hold on
        matHalfSize = floor(size(SegmentedChildrenObjectImage)/2);
        VisualizationSize = floor(size(SegmentedChildrenObjectImage)/6);
        
        subplot(2,2,1);
        %     matDistancesToOutline_Temp = [0;matDistancesToOutline];
        matDistancesToOutline_Temp = [0;matLocMeasurementsZScored(:,1)];
        ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
        [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
        ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
        
        
        ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
        imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
            matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [max(-3,quantile(ChildenColorCode_Temp(:),0.001)) max(1,quantile(ChildenColorCode_Temp(:),0.999))])
        
        matColorMap = jet; matColorMap(1,:) = 1;
        colormap(matColorMap)
        colorbar
        title('Min DistToCell Membr.')
        
        subplot(2,2,3);
        %     matDistancesToOutline_Temp = [0;matDisToOutlineThirdLine];
        matDistancesToOutline_Temp = [0;matLocMeasurementsZScored(:,5)];
        
        ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
        [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
        ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
        
        
        ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
        imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
            matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [max(-3,quantile(ChildenColorCode_Temp(:),0.001)) max(1,quantile(ChildenColorCode_Temp(:),0.999))])
        matColorMap = jet; matColorMap(1,:) = 1;
        colormap(matColorMap)
        colorbar
        title('Centroid Tangent DistToCell Membr.')
        
        
        
        
        subplot(2,2,2);
        
        %     matDistancesToOutline_Temp = [0;matDisPerX1(:,1)];
        matDistancesToOutline_Temp = [0;matLocMeasurementsZScored(:,9)];
        
        ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
        [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
        ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
        
        
        ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
        imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
            matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [max(-3,quantile(ChildenColorCode_Temp(:),0.001)) max(1,quantile(ChildenColorCode_Temp(:),0.999))])
        matColorMap = jet; matColorMap(1,:) = 1;
        colormap(matColorMap)
        colorbar
        title('Radious for first percentage')
        
        
        
        subplot(2,2,4);
        
        %     matDistancesToOutline_Temp = [0;matMeanDis];
        matDistancesToOutline_Temp = [0;matLocMeasurementsZScored(:,6)];
        
        ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
        [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
        ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
        
        ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
        imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
            matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [max(-3,quantile(ChildenColorCode_Temp(:),0.001)) max(1,quantile(ChildenColorCode_Temp(:),0.999))])
        matColorMap = jet; matColorMap(1,:) = 1;
        colormap(matColorMap)
        colorbar
        title('Mean distance to spots')
        
        
    end
    drawnow
catch
    
end



