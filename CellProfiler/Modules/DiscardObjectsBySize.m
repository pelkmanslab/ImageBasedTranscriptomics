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
%   Berend Snijder
%   Stephan Daetwyler
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
%defaultVAR02 = 2000
MinAreaSize = str2double(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What is the minimum allowed object area at the border?
%defaultVAR03 = 1000
MinAreaSizeBorder = str2double(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How large is the expected border?
%defaultVAR04 = 50
BorderSize = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

TargetName = ObjectName;

LabelMatrixImage = CPretrieveimage(handles,['Segmented' ObjectName],ModuleName,'MustBeGray','DontCheckScale');

matObjectSizes = regionprops(LabelMatrixImage,'Area');
matObjectSizes = cat(1,matObjectSizes.Area);

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow
  
%first get position of objects, make a list of those objects 

positionofobject=handles.Measurements.Nuclei.Location{handles.Current.SetBeingAnalyzed};
atBorder=[];
imagesizeforborder=size(LabelMatrixImage);
for i=1:length(positionofobject(:,1))
   if find(positionofobject(i,1)<BorderSize | positionofobject(i,2)<BorderSize | positionofobject(i,1)>imagesizeforborder(1,1)-BorderSize | positionofobject(i,2)> imagesizeforborder(1,2)-BorderSize)
       atBorder = [atBorder;i];
   end
end

%make a distinction between objects at the border and not at the border for
%size comparison
Filter = [];
for index = 1:length(positionofobject(:,1))
    if find(index ==atBorder)
        if matObjectSizes(index) < MinAreaSizeBorder
        Filter = [Filter;index];
        end
    else
        if matObjectSizes(index)<MinAreaSize
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

FinalLabelMatrixImage = LookUpColumn(FinalLabelMatrixImage+1);
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


% if ~strcmp(SaveOutlines,'Do not save')
%     try handles.Pipeline.(SaveOutlines) = LogicalOutlines;
%     catch
%         error(['The object outlines were not calculated by the ', ModuleName, ' module so these images were not saved to the handles structure. Image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
%     end
% end