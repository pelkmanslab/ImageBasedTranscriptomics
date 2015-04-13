function handles = DiscardSinglePixelObjects(handles)

% Help for the Discard Objects By Single Pixels :
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Eliminates objects if they are less than 2 pixels big
% *************************************************************************
%
%
% See also MeasureObjectAreaShape, MeasureObjectIntensity, MeasureTexture,
% MeasureCorrelation, CalculateRatios, and MeasureObjectNeighbors modules.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
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
%defaultVAR02 = 2
MinAreaSize = str2double(handles.Settings.VariableValues{CurrentModuleNum,2});

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

Filter = find(matObjectSizes < MinAreaSize);
FinalLabelMatrixImage = LabelMatrixImage;
for i=1:numel(Filter)
    FinalLabelMatrixImage(FinalLabelMatrixImage == Filter(i)) = 0;
end

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