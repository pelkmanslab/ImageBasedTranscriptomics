function handles = IdentifySpots2D(handles)
% Help for the IdentifySpots2D module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Detects spots as destribed by Battich et al., 2013.
% ***********************************************
% Will Determine Spots in 2D Image stacks after Laplacian Of Gaussian (LoG)
% enhancing of spots. Many of the input arguments are optional. Note that
% while an external script has to be run in order to choose robust values,
% manual selection of the parameters can often yield good estimates, if
% the signal is clear enough.
%
% WHAT DID YOU CALL THE IMAGES YOU WANT TO PROCESS?
% Object detection should be done on this image.
%
% HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED PRIOR TO DEBLENDING?
% This is the name of the the spots identified after thresholding the LoG
% image.
%
% HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED AFTER DEBLENDING?
% Optional. Deblending can be done after spot detection to separate close
% objects. The algorithm is based upon SourceExtractor. To skip this step,
% insert / as name of the object.
%
% OBJECTSIZE
% This value corresponds to the approximate size of you spots. It should
% be their diameter in pixels. The LoG will use a mask of this size to
% enhance radial signal of that size. Note that in practice the specific value
% does not affect the number of spots, if spots are bright (eg. pixel size 5
% or 6).
%
% INTENSITY QUANTA PER IMAGE
% Prior to spot detection the images are rescaled according to their
% intensity. Since the specific value of minimal and maximal intensities
% are frequently not robust across multiple images, intensity quantile are
% used instead. [0 1] would correspond to using the single dimmest pixel
% for minimal intensity and the single brightest pixel for maximal
% intensity. [0.01 0.90] would mean that the minimum intensity is derived
% from the pixel, which is the 1% brightest pixel of all and that the
% maximum intensity is derived from the pixel, which is the 90% brightest
% pixel .
%
% INTENSITY BORERS FOR INTENSITY RESCALING OF IMAGES
% Most extreme values that the image intensity minimum and image intensity
% maximum (as defined by the quanta) are allowed to have
% [LowestPossibleGreyscaleValueForImageMinimum
% HighestPossibleGreyscaleValueForImageMinimum
% LowestPossibleGreyscaleValueForImageMaximum
% HighestPossibleGreyscaleValueForImageMaximum]
% To ignore individual values, place a NaN.
% Note that these parameters very strongly depend upon the variability of
% your illumination source. When using a robust confocal microscope you can
% set the lowest and highest possible values to values,  which are very
% close (or even identical). If your light source is variable during the
% acquisition (which can be the case with Halogen lamps) you might choose
% less strict borders to detect spots of varying intensites.
%
% THRESHOLD OF SPOT DETECTION
% This is the threshold value for spot detection. The higher it is the more
% stringent your spot detection is. Use external script to determine a
% threshold where the spot number is robust against small variations in the
% threshold.
%
% HOW MANY STEPS OF DEBLENDING DO YOU WANT TO DO?
% The amount of deblending steps, which are done. The higher it is the less
% likely it is that two adjacent spots are not separated. The default of 30
% works very well (and we did not see improvement on our images with higher
% values). Note that the number of deblending steps is the main determinant
% of computational time for this module.
%
% WHAT IS THE MINIMAL INTENSITY OF A PIXEL WITHIN A SPOT?
% Minimal greyscale value of a pixel, which a pixel has to have in order to
% be recognized to be within a spot. Opitonal argument to make spot
% detection even more robust against very dim spots. In practice, we have
% never observed that this parameter would have any influence on the spot
% detection. However, you might include it as an additional safety measure.
%
% WHICH IMAGE DO YOU WANT TO USE AS A REFERENCE FOR SPOT BIAS CORRECTION?
% Here you can name a correction matrix which counteracts bias of the spot
% correction across the field of view. Note that such a correction matrix
% has to be loaded previously by a separate module, such as
% LOADSINGLEMATRIX
%
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
% The design of this module largely follows a IdentifyPrimLoG2 by
% Baris Sumengen.
%
% $Revision: 1889 $


drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How do you want to call the objects identified PRIOR to deblending?
%defaultVAR02 = PreSpots
%infotypeVAR02 = objectgroup indep
iPreObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = How do you want to call the objects identified AFTER deblending?
%defaultVAR03 = /
%infotypeVAR03 = objectgroup indep
iPostObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = ObjectSize
%defaultVAR04 = 6
iHsize = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Intensity Quanta Per Image
%defaultVAR05 = [0.01 0.99]
iImgLimes = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Intensity borders for intensity rescaling of images
%[MinOfMinintens MaxOfMinintens MinOfMaxintens MaxOfMaxintens]
%defaultVAR06 = [NaN 120 500 NaN]
iRescaleThr = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Threshold of Spot Detection
%defaultVAR07 = 0.01
iDetectionThr = char(handles.Settings.VariableValues{CurrentModuleNum,7});


%textVAR08 = How many Steps of Deblending do you want to do?
%defaultVAR08 = 0
iDeblendSteps = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = What is the minimal intensity of a pixel within a spot?
%defaultVAR09 = /
iObjIntensityThr = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = Do you want to perform spot bias correction?
%choiceVAR10 = No
%choiceVAR10 = Yes
iDoBiasCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = Which image do you want to use as a reference for spot bias correction?
%infotypeVAR11 = imagegroup
iCorrectionName = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%%VariableRevisionNumber = 5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHECK INPUT   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter Size
try
    iHsize = str2double(iHsize);
catch errFilterSize
    error(['Image processing was canceled in the ', ModuleName, ' module because the object size could not be converted to a number.'])
end

if iHsize<=2
    error(['Image processing was canceled in the ', ModuleName, ' module because the object size was too small. Has to be at least 3'])
end

% Intensity Quanta Of Image
[isSafe iImgLimes]= inputVectorsForEvalCP3D(iImgLimes,true);
if isSafe ==false
    error(['Image processing was canceled in the ', ModuleName, ' module because Intensity Quanta per Image contain forbidden characters.'])
end

% Rescale Thresholds
[isSafe iRescaleThr]= inputVectorsForEvalCP3D(iRescaleThr,true);
if isSafe ==false
    error(['Image processing was canceled in the ', ModuleName, ' module because Rescaling Boundaries contain forbidden characters.'])
end

% Detection Threshold
try
    iDetectionThr = str2double(iDetectionThr);
catch errDetectionThr
    error(['Image processing was canceled in the ', ModuleName, ' module because the Detection Threshold could not be converted to a number.'])
end

% Deblend Threshold
try
    iDeblendSteps = str2double(iDeblendSteps);
catch errDeblendDetection
    error(['Image processing was canceled in the ', ModuleName, ' module because the Stepsize for deblending could not be converted to a number.'])
end


if iPostObjectName == '/'
    if iDeblendSteps > 0
        error(['Image processing was canceled in the ', ModuleName, ' module because no Name for the Objects after deblending were defined.'])
    end
else
    if iDeblendSteps <0
        error(['Image processing was canceled in the ', ModuleName, ' module because the amount of deblending steps were not defined.'])
    end
end

% Bias correction
if strcmp(iDoBiasCorrection,'Yes')
    bnDoBiasCorrection = true;
else
    bnDoBiasCorrection = false;
end


% Initiate Bias Correction
if bnDoBiasCorrection == true
    DetectionBias =  handles.Pipeline.(iCorrectionName);
else
    DetectionBias = [];
end

% Initiate Settings for deblending
Options.ObSize = iHsize;
Options.limQuant = eval(iImgLimes);
Options.RescaleThr = eval(iRescaleThr);
Options.ObjIntensityThr = [];
Options.closeHoles = false;
Options.ObjSizeThr = [];
Options.ObjThr = iDetectionThr;
Options.StepNumber = iDeblendSteps;
Options.numRatio = 0.20;
Options.doLog = 0;
Options.DetectBias = DetectionBias;

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
Image = double(CPretrieveimage(handles,ImageName,ModuleName,'DontCheckColor','DontCheckScale')).*65535;     % convert to scale used for spotdetection
op = fspecialCP3D('2D LoG',iHsize);         % force 2D filter

% Object intensity Threshold
if iObjIntensityThr == '/'
    iObjIntensityThr = [];
elseif strcmp(iObjIntensityThr,'auto')
    iObjIntensityThr = multithresh(Image).*0.7;% be sensitive enough!
else
    try
        iObjIntensityThr = str2double(iObjIntensityThr);
    catch errObjIntensityThr
        error(['Image processing was canceled in the ', ModuleName, ' module because the Stepsize for deblending could not be converted to a number.'])
    end
end

% Detect objects, note that input vectors are eval'ed
[ObjCount{1} SegmentationCC{1} FiltImage] = ObjByFilter(Image,op,iDetectionThr,eval(iImgLimes),eval(iRescaleThr),iObjIntensityThr,true,[],DetectionBias);
% Convert to CP1 standard: labelmatrix
MatrixLabel{1} = double(labelmatrix(SegmentationCC{1}));

% Security check, if conversion is correct
if max(MatrixLabel{1}(:)) ~= ObjCount{1}
    error(['Image processing was canceled in the ', ModuleName, ' module because conversion of format of segmentation was wrong. Contact Thomas.'])
end

% Deblend objects
if iDeblendSteps > 0        % Only do deblending, if number of iterations was defined
    MatrixLabel{2} = SourceExtractorDeblend(Image,SegmentationCC{1},FiltImage,Options);
    ObjCount{2} = max(MatrixLabel{2}(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MEASUREMENTS and SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectName{1} = iPreObjectName;
if iDeblendSteps > 0
    numObjects = 2;
    ObjectName{2} = iPostObjectName;
else
    numObjects = 1;
end

for o = 1:numObjects % Cycle through object groups that should be saved (pre and post deblending)
    
    %%% Save Segmentation to Pipeline
    fieldname = ['Segmented', ObjectName{o}];
    handles.Pipeline.(fieldname) = MatrixLabel{o};
    
    fieldname = ['SmallRemovedSegmented', ObjectName{o}];
    handles.Pipeline.(fieldname) = MatrixLabel{o};
    
    fieldname = ['UneditedSegmented', ObjectName{o}];
    handles.Pipeline.(fieldname) = MatrixLabel{o};
    
    %%% Saves the ObjectCount, i.e. the number of segmented objects:
    if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
        handles.Measurements.Image.ObjectCountFeatures = {};
        handles.Measurements.Image.ObjectCount = {};
    end
    column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName{o}));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName{o}};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount{o};
    
    % Save Centroid
    handles.Measurements.(ObjectName{o}).LocationFeatures = {'CenterX','CenterY'};
    
    Centroid = [0 0];
    if ObjCount{o} ~= 0 % determine centroid, if at least one object
        tmp = regionprops(MatrixLabel{o},'Centroid');
        Centroid = cat(1,tmp.Centroid);
    end
    handles.Measurements.(ObjectName{o}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
end

%%%%%%%%%%%%%%%%%%%
%%% DISPLAY %%%%%%%
%%%%%%%%%%%%%%%%%%%

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless == false
        % Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        
        subplot(2,1,1)      % Subplot with input image
        CPimagesc(Image,handles);
        title([ImageName, ' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        
        
        subplot(2,1,2)
        switch numObjects
            case 1
                bwImage = MatrixLabel{1}>0;
            case 2  % in case that deblending was done, dilate by 2 pixels to help visualization
                bwImage = imdilate(MatrixLabel{2}>0, strel('disk', 2));
        end
        
        r = (Image - min(Image(:))) / quantile(Image(:),0.995);
        g = (Image - min(Image(:))) / quantile(Image(:),0.995);
        b = (Image - min(Image(:))) / quantile(Image(:),0.995);
        
        r(bwImage) = max(r(:));
        g(bwImage) = 0;
        b(bwImage) = 0;
        visRGB = cat(3, r, g, b);
        f = visRGB <0;
        visRGB(f)=0;
        f = visRGB >1;
        visRGB(f)=1;
        
        
        CPimagesc(visRGB, handles);
        switch numObjects
            case 1
                title([ObjectName{1} ' (no deblending) Total count' num2str(ObjCount{1})]);
            case 2
                title([ObjectName{2} ' Total count' num2str(ObjCount{2}) ' (after deblending) ' num2str(ObjCount{1}) ' (before deblending)']);
        end
    end
end



end
