function handles = IdentifyPrimaryIterative(handles)

% Help for IdentifyPrimaryIterative
% Category: Object Processing
%
%
% DESCRIPTION:
% Primary identification of objects (nuclei) based on intensity thresholding 
% and subsequent separation of clumped objects along watershed lines between
% concave regions.
%
% DETAILS:
% Pixels belonging to nuclei objects can be easily separated form background 
% pixels by thresholding an image of a nuclei-specific stain such as DAPI.
% However, this often results in clumps of several individual objects, 
% because a single, image-wide threshold value is generally not
% sufficient to nicely separate objects, which lie very close to each other.
% Such clumped objects have a distinct morphology. Compared to individual objects, 
% which are more or less round and lie within a certain size range,
% clumped objects are relatively large and display multiple concave regions. 
% The intersection of individual objects is most likely a line connecting two 
% concave regions. This separating cut line can be found using the watershed 
% algorithm. By restricting watershed lines to the area between two concave regions, 
% very stable segmentation results can be achieved.
% The module processes the input image as follows:
% 1) Initial objects are identified by simple thresholding.
% 2) Clumped objects are selected on the basis of size and shape features:
%    area, solidity, and form factor.
% 3) The perimeter of selected objects is analyzed and concave 
%    regions along the boundary of objects are identified.
% 4) Watershed lines connecting two concave regions are determined.
% 5) All possible cuts along the selected watershed lines are considered and features
%    of each cut line (intensity along the line, angle between concave regions)
%    as well as features of the resulting objects (area/shape) are measured.
% 6) An "optimal" cut line is finally chosen by minimizing a cost function that
%    evaluates the measured features such that the resulting objects have a minimal
%    size and are as round as possible, while the separating line is as straight
%    and short as possible and the intensity along the line as low as possible. 
% Note that once an object has been selected for cutting and concave regions
% have been identified, a cut is inevitably made! 
% You can control the selection of clumped objects and the identification of 
% concave regions by setting the corresponding parameters (see below).
% Test modes are available for both steps that allow choosing parameter values 
% in a visually assisted way.
% 
%
% PARAMETERS:
% Object name:
% The name you would like to give the objects identified by this module.
%
% Image name:
% The name of the input image in which primary objects should be identified.
% 
% Threshold correction factor:
% When the threshold is calculated automatically, it may consistently be
% too stringent or too lenient. You may need to enter an adjustment factor,
% which you empirically determine suitable for your images. The number 1
% means no adjustment, 0 to 1 makes the threshold more lenient, and greater
% than 1 (e.g. 1.3) makes the threshold more stringent. The thresholding method
% used by this module (Otsu algorithm) inherently assumes that 50% of the image
% is covered by objects. If a larger percentage of the image is covered, the
% method will give a slightly biased threshold that may have to be
% corrected using a threshold correction factor.
%
% Lower and upper bounds on threshold:
% Can be used as a safety precaution when the threshold is calculated
% automatically. For example, if there are no objects in the field of view,
% the automatic threshold will be unreasonably low. In such cases, the
% lower bound you enter here will override the automatic threshold.
% 
% Cutting passes: 
% Each pass, only one cut per concave region is allowed, possibly making it
% necessary to perform additional cutting passes to separate clumps of 
% more than two objects.
%
% Debug mode:
% Separation of clumped objects is done on small sub-images.
% By activating the debug mode, you can visually follow steps 2-5 of the algorithm
% outlined above and check whether your parameter settings have the desired effect,
% i.e. whether the correct regions and lines are selected for each object.
%
% Object selection:
% Limits for solidity, form factor, and upper and lower size of objects to be 
% selected for cutting. Determine optimal values via test mode.
% 
% Test mode for object selection: 
% Displays solidity, area, and (transformed) form factor values for each primary
% object identified by thresholding. Pick values from images to fine tune settings.
%
% Perimeter analysis: 
% Parameters for detection of concave regions. Determine optimal value via test mode.
% Window size:
% Sliding window for calculating the curvature of objects. Large values result
% in more continuous, smoother but maybe less precise regions, while small values
% give more precise, but smaller and less continuous regions.
% Max equivalent radius: 
% Maximum equivalent radius of a concave region to be eligible for cutting.
% Higher values increase sensitivity and result in more cut options.
% Min equivalent segment: 
% Minimum equivalent circular fragment (degree) of a concave region to be 
% eligible for cutting. Lower values increase sensitivity and result in more cut options.
%
% Test mode for perimeter analysis: 
% Displays curvature, convex/concave, equivalent radius and segment of each object.
% Pick values from images to fine tune settings.
%
% DEPENDENCIES:
% PerimeterAnalysis.m
% PerimeterWatershedSegmentation.m
%
% AUTHORS:
%  Markus Herrmann
%  Nicolas Battich
%  Thomas Stoeger
%  Anatol Schwab
%
% (c) Pelkmans Lab 2015
%
% $Revision: 1879 $



%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the objects identified by this module?
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the intensity image that should be used for object identification?
%infotypeVAR02 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Intensity thresholding: Threshold correction factor
%defaultVAR03 = 1
ThresholdCorrection = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,3}));

%textVAR04 = Intensity thresholding: Lower and upper bounds on threshold, in the range [0,1]
%defaultVAR04 = 0,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Cutting passes (0 = no cutting)
%defaultVAR05 = 2
CuttingPasses = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Debug mode: Show individual steps while iterating over clumped objects
%choiceVAR06 = Off
%choiceVAR06 = On
DebugMode = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 = Object selection: Maximal SOLIDITY of objects, which should be cut (1 = solidity independent)
%defaultVAR07 = 0.92
SolidityThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = Object selection: Minimal FORM FACTOR of objects, which should be cut (0 = form factor independent)
%defaultVAR08 = 0.40
FormFactorThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = Object selection: Maximal AREA of objects, which should be cut (0 = area independent)
%defaultVAR09 = 50000
UpperSizeThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = Object selection: Minimal AREA of objects, which should be cut (0 = area independent)
%defaultVAR10 = 5000
LowerSizeThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Minimal area that cut objects should have.
%defaultVAR11 = 2000
LowerSizeCutThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,11}));

%textVAR12 = Test mode for object selection: solidity, area, form factor
%choiceVAR12 = No
%choiceVAR12 = Yes
TestMode2 = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = Perimeter analysis: SLIDING WINDOW size for curvature calculation
%defaultVAR13 = 8
WindowSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,13}));

%textVAR14 = Perimeter analysis: FILTER SIZE for smoothing objects
%defaultVAR14 = 1
smoothingDiskSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,14}));

%textVAR15 = Perimeter analysis: Maximum concave region equivalent RADIUS
%defaultVAR15 = 30
PerimSegEqRadius = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,15}));

%textVAR16 = Perimeter analysis: Minimum concave region equivalent CIRCULAR SEGMENT (degree)
%defaultVAR16 = 6
PerimSegEqSegment = degtorad(str2double(char(handles.Settings.VariableValues{CurrentModuleNum,16})));

%textVAR17 = Test mode for perimeter analysis: overlay curvature etc. on objects
%choiceVAR17 = No
%choiceVAR17 = Yes
TestMode = char(handles.Settings.VariableValues{CurrentModuleNum,17});
%inputtypeVAR17 = popupmenu

%%%VariableRevisionNumber = 15



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IMAGES FROM HANDLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigImage = handles.Pipeline.(ImageName);


%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%

Threshold = 'Otsu Global';
pObject = '10%';

%%% Chooses the first word of the method name (removing 'Global' or 'Adaptive').
ThresholdMethod = strtok(Threshold);
%%% Checks if a custom entry was selected for Threshold, which means we are using an incoming binary image rather than calculating a threshold.
if isempty(strmatch(ThresholdMethod,{'Otsu','MoG','Background','RobustBackground','RidlerCalvard','All','Set'},'exact'))
    %if ~(strncmp(Threshold,'Otsu',4) || strncmp(Threshold,'MoG',3) || strfind(Threshold,'Background') ||strncmp(Threshold,'RidlerCalvard',13) || strcmp(Threshold,'All') || strcmp(Threshold,'Set interactively'))
    if isnan(str2double(Threshold))
        GetThreshold = 0;
        BinaryInputImage = CPretrieveimage(handles,Threshold,ModuleName,'MustBeGray','CheckScale');
    else
        GetThreshold = 1;
    end
else
    GetThreshold = 1;
end

%%% Checks that the Min and Max threshold bounds have valid values
index = strfind(ThresholdRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
end
MinimumThreshold = ThresholdRange(1:index-1);
MaximumThreshold = ThresholdRange(index+1:end);


if GetThreshold
    OrigImage(OrigImage > quantile(OrigImage(:), 0.999)) = quantile(OrigImage(:), 0.999);
    [handles,OrigThreshold] = CPthreshold(handles,Threshold,pObject,MinimumThreshold,MaximumThreshold,ThresholdCorrection,OrigImage,ImageName,ModuleName,ObjectName);
else
    OrigThreshold = 0;
end

%%% Threshold intensity image
ThreshImage = zeros(size(OrigImage), 'double');
ThreshImage(OrigImage > OrigThreshold) = 1;

%%% Fill holes in objects
imInputObjects = imfill(double(ThreshImage),'holes');



if ~isempty(imInputObjects)
    
    %-------------------------------------------
    % Select objects in input image for cutting
    %-------------------------------------------
    
    imObjects = zeros([size(imInputObjects),CuttingPasses]);
    imSelected = zeros([size(imInputObjects),CuttingPasses]);
    imCutMask = zeros([size(imInputObjects),CuttingPasses]);
    imCut = zeros([size(imInputObjects),CuttingPasses]);
    imNotCut = zeros([size(imInputObjects),CuttingPasses]);
    objFormFactor = cell(CuttingPasses,1);
    objSolidity = cell(CuttingPasses,1);
    objArea = cell(CuttingPasses,1);
    cellPerimeterProps = cell(CuttingPasses,1);
    
    for i = 1:CuttingPasses
        
        if i==1
            imObjects(:,:,i) = imInputObjects;
        else
            imObjects(:,:,i) = imCut(:,:,i-1);
        end
        
        % Measure basic area/shape features
        props = regionprops(logical(imObjects(:,:,i)),'Area','Solidity','Perimeter');

        % Features used for object selection
        objSolidity{i} = cat(1,props.Solidity);
        objArea{i} = cat(1,props.Area);
        tmp = log((4*pi*cat(1,props.Area)) ./ ((cat(1,props.Perimeter)+1).^2))*(-1);%make values positive for easier interpretation of parameter values
        tmp(tmp<0) = 0;
        objFormFactor{i} = tmp;

        % Select objects based on these features (user defined thresholds)
        obj2cut = objSolidity{i} < SolidityThres & objFormFactor{i} > FormFactorThres ...
            & objArea{i} > LowerSizeThres & objArea{i} < UpperSizeThres;
        objNot2cut = ~obj2cut;
                    
        objSelected = zeros(size(obj2cut));
        objSelected(obj2cut) = 1;
        objSelected(objNot2cut) = 2;
        imSelected(:,:,i) = rplabel(logical(imObjects(:,:,i)),[],objSelected);
        
        % Create mask image with objects selected for cutting
        imObj2Cut = zeros(size(OrigImage));
        imObj2Cut(imSelected(:,:,i)==1) = 1;
        
        % Store remaining objects that are omitted from cutting
        tmp = zeros(size(OrigImage));
        tmp(imSelected(:,:,i)==2) = 1;
        imNotCut(:,:,i) = logical(tmp);
        
        
        %-------------
        % Cut objects
        %-------------
        
        % Smooth image
        SmoothDisk = getnhood(strel('disk',smoothingDiskSize,0));%minimum that has to be done to avoid problems with bwtraceboundary
        imObj2Cut = bwlabel(imdilate(imerode(imObj2Cut,SmoothDisk),SmoothDisk));
        
        % Separate clumped objects along watershed lines

        % Note: PerimeterAnalysis cannot handle holes in objects (we may
        % want to implement this in case of big clumps of many objects).
        % Sliding window size is linked to object size. Small object sizes
        % (e.g. in case of images acquired with low magnification) limits
        % maximal size of the sliding window and thus sensitivity of the
        % perimeter analysis.
        
        % Perform perimeter analysis
        cellPerimeterProps{i} = PerimeterAnalysis(imObj2Cut,WindowSize);

        % This parameter limits the number of allowed concave regions.
        % It can serve as a safety measure to prevent runtime problems for
        % very complex objects.
        % This could become an input argument in the future!?
        numRegionTheshold = 30;
        
        % Perform the actual segmentation
        if strcmp(DebugMode, 'On')
            imCutMask(:,:,i) = PerimeterWatershedSegmentation(imObj2Cut,OrigImage,cellPerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,LowerSizeCutThres, numRegionTheshold, 'debugON');
        else
            imCutMask(:,:,i) = PerimeterWatershedSegmentation(imObj2Cut,OrigImage,cellPerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,LowerSizeCutThres, numRegionTheshold);
        end
        imCut(:,:,i) = bwlabel(imObj2Cut.*~imCutMask(:,:,i));
        
        
        %------------------------------
        % Display intermediate results
        %------------------------------
        
        drawnow
        
        % Create overlay images
        imOutlineShapeSeparatedOverlay = OrigImage;
        B = bwboundaries(imCut(:,:,i));
        imCutShapeObjectsLabel = label2rgb(bwlabel(imCut(:,:,i)),'jet','k','shuffle');
        
        % GUI
        tmpSelected = (imSelected(:,:,i));
        ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        subplot(2,2,2), CPimagesc(logical(tmpSelected==1),handles),
        title(['Cut lines on selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        hold on
        L = bwboundaries(logical(imCutMask(:,:,i)), 'noholes');
        for l = 1:length(L)
            line = L{l};
            plot(line(:,2), line(:,1), 'r', 'LineWidth', 3);
        end
        hold off
        freezeColors
        subplot(2,2,1), CPimagesc(imSelected(:,:,i),handles), colormap('jet'),
        title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        freezeColors
        subplot(2,2,3), CPimagesc(imOutlineShapeSeparatedOverlay,handles),
        title(['Outlines of Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
        end
        hold off
        freezeColors
        subplot(2,2,4), CPimagesc(imCutShapeObjectsLabel,handles),
        title(['Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        freezeColors
        
    end
    
    %-----------------------------------------------
    % Combine objects from different cutting passes
    %-----------------------------------------------
    
    imCut = logical(imCut(:,:,CuttingPasses));
    
    if ~isempty(imCut)
        imErodeMask = bwmorph(imCut,'shrink',inf);
        imDilatedMask = IdentifySecPropagateSubfunction(double(imErodeMask),OrigImage,imCut,1);
    end
    
    imNotCut = logical(sum(imNotCut,3));% Retrieve objects that were not cut
    imFinalObjects = bwlabel(logical(imDilatedMask + imNotCut));
    
else
    
    cellPerimeterProps = {};
    imFinalObjects = zeros(size(imInputObjects));
    imObjects = zeros([size(imInputObjects),CuttingPasses]);
    imSelected = zeros([size(imInputObjects),CuttingPasses]);
    imCutMask = zeros([size(imInputObjects),CuttingPasses]);
    imCut = zeros([size(imInputObjects),CuttingPasses]);
    imNotCut = zeros([size(imInputObjects),CuttingPasses]);
    objFormFactor = cell(CuttingPasses,1);
    objSolidity = cell(CuttingPasses,1);
    objArea = cell(CuttingPasses,1);
    cellPerimeterProps = cell(CuttingPasses,1);
    
end


%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%
%%%%%%%%%%%%%%%%%%%%%

drawnow

% Create overlay images
imOutlineShapeSeparatedOverlay = OrigImage;
B = bwboundaries(logical(imFinalObjects),'holes');
imCutShapeObjectsLabel = label2rgb(bwlabel(imFinalObjects),'jet','k','shuffle');

% GUI
% CPfigure(handles,'Image',ThisModuleFigureNumber);
imDisplay = imSelected(:,:,1);
subplot(2,2,2), CPimagesc(logical(imDisplay),handles),
title(['Cut lines on selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
hold on
L = bwboundaries(logical(sum(imCutMask, 3)), 'noholes');
for i = 1:length(L)
    line = L{i};
    plot(line(:,2), line(:,1), 'r', 'LineWidth', 3);
end
hold off
freezeColors
subplot(2,2,1), CPimagesc(imDisplay,handles), colormap('jet'),
title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
freezeColors
subplot(2,2,3), CPimagesc(imOutlineShapeSeparatedOverlay,handles),
title(['Outlines of Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
end
hold off
freezeColors
subplot(2,2,4), CPimagesc(imCutShapeObjectsLabel,handles),
title(['Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
freezeColors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmpi(savePDF, 'Yes')
%     % Save figure als PDF
%     strFigureName = sprintf('%s_Iteration%d',mfilename,handles.Current.SetBeingAnalyzed);
%     strPDFdir = strrep(handles.Current.DefaultOutputDirectory, 'BATCH', 'PDF');
%     if ~fileattrib(strPDFdir)
%         mkdir(strPDFdir);
%         fprintf('%s: creating directory: ''%s''\n', mfilename, strPDFdir);
%     end
%     gcf2pdf(strPDFdir, strFigureName)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Plot shape analysis data
if strcmp(TestMode,'Yes')
    if ~isempty(cellPerimeterProps)
        for h = 1:CuttingPasses
            imCurvature = zeros(size(OrigImage),'double');
            imConvexConcave = zeros(size(OrigImage),'double');
            imAngle = zeros(size(OrigImage),'double');
            imRadius = zeros(size(OrigImage),'double');
            for i = 1:length(cellPerimeterProps{h})
                matCurrentObjectProps = cellPerimeterProps{h}{i};%get current object
                imConcaveRegions = bwlabel(matCurrentObjectProps(:,11)==-1);
                imConvexRegions = bwlabel(matCurrentObjectProps(:,11)==1);
                AllRegions = imConcaveRegions+(max(imConcaveRegions)+imConvexRegions).*(imConvexRegions>0);%bwlabel only works binary, therefore label convex, concave seperately, then merger labels
                NumRegions = length(setdiff(unique(AllRegions),0));
                for j = 1:size(matCurrentObjectProps,1)%loop over all pixels of object to plot general properties
                    imCurvature(matCurrentObjectProps(j,1),matCurrentObjectProps(j,2)) = matCurrentObjectProps(j,9);
                    imConvexConcave(matCurrentObjectProps(j,1),matCurrentObjectProps(j,2)) = matCurrentObjectProps(j,11);
                end
                for k = 1:NumRegions%loop over all regions to plot region specific properties
                    matCurrentRegionProps = matCurrentObjectProps(AllRegions==k,:);%get current region
                    NormCurvature = matCurrentRegionProps(:,9);
                    CurrentEqAngle = sum(NormCurvature);
                    CurrentEqRadius = length(NormCurvature)/sum(NormCurvature);
                    for L = 1:size(matCurrentRegionProps,1)%loop over all pixels in region
                        imRadius(matCurrentRegionProps(L,1),matCurrentRegionProps(L,2)) = CurrentEqRadius;
                        imAngle(matCurrentRegionProps(L,1),matCurrentRegionProps(L,2)) = radtodeg(CurrentEqAngle);
                    end
                end
            end
            CPfigure('Tag',strcat('ShapeAnalysisPass',num2str(h)));
            
            subplot(2,2,1);
            CPimagesc(imCurvature,handles);
            title(['Curvature image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
            
            subplot(2,2,2);
            %problem with the CP image range scaling hack: while CPimagesc would
            %accept the range as an argument, 'Open in new window' will ignore
            %it. therefore the function has to be tricked somehow! solution:
            %make rgb image with each channel binary
            RGBConvexConcaveImage = cat(3,(imConvexConcave==1),(imConvexConcave==-1),zeros(size(imConvexConcave)));
            CPimagesc(RGBConvexConcaveImage,handles);
            title(['Convex concave image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
            
            subplot(2,2,3);
            CPimagesc(imAngle,handles);
            title(['Equivalent angle (degree) image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
            
            subplot(2,2,4);
            CPimagesc(imRadius,handles);
            title(['Equivalent radius, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
        end
    end
end

% Plot area/shape feature data
if strcmp(TestMode2,'Yes')
    if ~isempty(cellPerimeterProps)
        for h = 1:CuttingPasses
            imSolidity = rplabel(logical(imObjects(:,:,h)), [], objSolidity{h});
            imFormFactor = rplabel(logical(imObjects(:,:,h)), [], objFormFactor{h});
            imArea = rplabel(logical(imObjects(:,:,h)), [], objArea{h});
            
            % could be nicely done with cbrewer() but stupid 'freezeColors'
            % erases the indices!!! note that colorbars could be preserved
            % with 'cbfreeze'
            %             cmapR = cbrewer('seq', 'Reds', 9);
            %             cmapG = cbrewer('seq', 'Greens', 9);
            %             cmapB = cbrewer('seq', 'Blues', 9);
            
            CPfigure('Tag','Features for object selection')
            subplot(2,2,1), CPimagesc(imSolidity,handles), %colormap(cmapR),
            title(['Solidity of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            subplot(2,2,2), CPimagesc(imFormFactor,handles), %colormap(cmapG),
            title(['Form factor of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            subplot(2,2,3), CPimagesc(imArea,handles), %colormap(cmapB),
            title(['Area of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            subplot(2,2,4), CPimagesc(imSelected(:,:,h),handles), colormap('jet'),
            title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE DATA TO HANDLES STRUCTURE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldname = ['UneditedSegmented',ObjectName];%not edited for size or edge
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['SmallRemovedSegmented',ObjectName];%for IdentifySecondary.m
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['Segmented',ObjectName];%final label image
handles.Pipeline.(fieldname) = imFinalObjects;

%%% Saves location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(imFinalObjects,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves ObjectCount, i.e. number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(imFinalObjects(:));


end

