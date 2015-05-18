function handles = DiscardObjectsBySize(handles)

% Help for the Discard Objects by Size module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Eliminates objects if they are less than a given size, and treats the
% border of the images differently
% *************************************************************************
%
% This module calculates the size of objects and discard objects below a
% critical size defined by the user. It also gives the probability to keep
% larger objects at the border.
%
% Authors:
%   Stephan Daetwyler
%   Thomas Stoeger
%   Berend Snijder
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What objects should be discarded?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What is the minimum allowed object area?
%defaultVAR02 = 900
MinAreaSize = str2double(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Do you want to set a different minimum allowed object area for objects touching the border?
%choiceVAR03 = No
%choiceVAR03 = Yes
useAdditionalThresholdAtBorder = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = If yes, enter the number of pixels that define the border region.
%defaultVAR04 = 50
BorderSize = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = If yes, what is the minimum allowed object area for objects within the border region?
%defaultVAR05 = 0
MinAreaSizeBorder = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

TargetName = ObjectName;

LabelMatrixImage = CPretrieveimage(handles,['Segmented' ObjectName],ModuleName,'MustBeGray','DontCheckScale');
UseAsLabelInCaseNoBackgroundPresent = LabelMatrixImage;

if any(LabelMatrixImage(:) == 0)
    originalSegmentationHasBackground = true;
else
    originalSegmentationHasBackground = false;
end

matObjectSizes = regionprops(LabelMatrixImage,'Area');
matObjectSizes = cat(1,matObjectSizes.Area);

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if ~any(LabelMatrixImage(:) > 0)
    FinalLabelMatrixImage = LabelMatrixImage;
else % only process segmentation, if at least one object is present
    
    
    % Input check
    if isequal(useAdditionalThresholdAtBorder,'No')
        useAdditionalThresholdAtBorder = false;
    elseif isequal(useAdditionalThresholdAtBorder,'Yes')
        useAdditionalThresholdAtBorder = true;
    else
        error('Only options, which are allowed for useAdditionalThresholdAtBorder are No and Yes!')
    end
    
    % If user does not want a separtate size filter at border (now: default; in
    % past: not default), set minimal area at border to same value as minimal
    % area away from from border region (and thus center of image site)
    if useAdditionalThresholdAtBorder == false
        MinAreaSizeBorder = MinAreaSize;
        BorderSize = 50;
    end
    
    % Determine object centroids that are close to the border of the image
    positionofobject=handles.Measurements.Nuclei.Location{handles.Current.SetBeingAnalyzed};
    
    atBorder=[];
    imagesizeforborder=size(LabelMatrixImage);
    for i=1:length(positionofobject(:,1))
        if find(positionofobject(i,1)<BorderSize | positionofobject(i,2)<BorderSize | positionofobject(i,1)>imagesizeforborder(1,1)-BorderSize | positionofobject(i,2)> imagesizeforborder(1,2)-BorderSize)
            atBorder = [atBorder;i];
        end
    end
    
    % Identify objects, whose area is below the requested minimal area
    Filter = [];
    for index = 1:length(positionofobject(:,1))
        if find(index ==atBorder)
            if matObjectSizes(index) < MinAreaSizeBorder
                Filter = [Filter;index];
            end
        else
            if matObjectSizes(index) < MinAreaSize
                Filter = [Filter;index];
            end
        end
    end
    clear atBorder
    clear positionofobject
    
    FinalLabelMatrixImage = LabelMatrixImage;
    for i=1:numel(Filter)
        FinalLabelMatrixImage(FinalLabelMatrixImage == Filter(i)) = 0;
    end
    
    clear Filter
    
    x = sortrows(unique([LabelMatrixImage(:) FinalLabelMatrixImage(:)],'rows'),1);
    x(x(:,2)>0,2)=1:sum(x(:,2)>0);
    LookUpColumn = x(:,2);
    
    if originalSegmentationHasBackground == true % default / CP's original code
        FinalLabelMatrixImage = LookUpColumn(FinalLabelMatrixImage+1);
    else
        % [TS 150504: There is a rare bug (also in CP's original
        % DiscardSinglePixelObjects) that causes a crash if no pixel belongs to
        % a backround; I have tested the seemingly obvious fix to remove the +1
        % so that the upper command is
        % FinalLabelMatrixImage = LookUpColumn(FinalLabelMatrixImage);
        % but unfortunately noted that this can cause other bugs, if a small
        % object is completely embedded in one large object that fills all
        % remaining pixels of the site
        FinalLabelMatrixImage = UseAsLabelInCaseNoBackgroundPresent;
        warning('%s: Preserve original segmentation since every pixel belongs to an object.\n', mfilename);
    end
    
end

%
% %%% Note: these outlines are not perfectly accurate; for some reason it
% %%% produces more objects than in the original image.  But it is OK for
% %%% display purposes.
% %%% Maximum filters the image with a 3x3 neighborhood.
% MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
% %%% Determines the outlines.
% IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
% %%% Converts to logical.
% warning off MATLAB:conversionToLogical
% LogicalOutlines = logical(IntensityOutlines);
% warning on MATLAB:conversionToLogical
% %%% Determines the grayscale intensity to use for the cell outlines.
% LineIntensity = max(OrigImage(:));
% %%% Overlays the outlines on the original image.
% ObjectOutlinesOnOrigImage = OrigImage;
% ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless == false
        
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        
        %%% A subplot of the figure window is set to display the original
        %%% image.
        subplot(2,2,1);
        ColoredLabelMatrixImage1 = CPlabel2rgb(handles,LabelMatrixImage);
        CPimagesc(ColoredLabelMatrixImage1,handles);
        title(sprintf('Input label matrix (%s): max = %d',ObjectName,max(LabelMatrixImage(:))));
        %%% A subplot of the figure window is set to display the label
        %%% matrix image.
        subplot(2,2,3);
        ColoredLabelMatrixImage2 = CPlabel2rgb(handles,FinalLabelMatrixImage);
        CPimagesc(ColoredLabelMatrixImage2,handles);
        title(sprintf('Output label matrix (%s): max = %d (filtering objects smaller than %d pixels)',ObjectName,max(FinalLabelMatrixImage(:)),MinAreaSize));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

handles.Pipeline.(['Segmented' TargetName]) = FinalLabelMatrixImage;

fieldname = ['SmallRemovedSegmented', ObjectName];
%%% Checks whether the image exists in the handles structure.
if isfield(handles.Pipeline, fieldname)
    handles.Pipeline.(['SmallRemovedSegmented' TargetName]) = handles.Pipeline.(['SmallRemovedSegmented',ObjectName]);
end

fieldname = ['UneditedSegmented',ObjectName];
%%% Checks whether the image exists in the handles structure.
if isfield(handles.Pipeline, fieldname)
    handles.Pipeline.(['UneditedSegmented' TargetName]) = handles.Pipeline.(['UneditedSegmented',ObjectName]);
end

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,TargetName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {['ObjectCount ' TargetName]};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(FinalLabelMatrixImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(TargetName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(FinalLabelMatrixImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];
end
handles.Measurements.(TargetName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves how many spots were discarded
handles.Measurements.(TargetName).DiscardedTooSmallObjectsFeatures = {'discarded'};
handles.Measurements.(TargetName).DiscardedTooSmallObjects(handles.Current.SetBeingAnalyzed) = {(max(LabelMatrixImage(:)) -  max(FinalLabelMatrixImage(:)))};
