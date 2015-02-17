function H = applyNBBSIllumCorrCP3D(Image,varargin)
%applyNBBSIllumCorrCP3D applies illumination according to NBBS method on 2D
%or 3D images. Correction can be done with a single correction reference
%for all z-planes of IMAGE, or with individual correction references for
%each z-plane
%
%   H = APPLYNBBSILLUMCORRCP3D(Image,STRINGARRAY) . STRINGARRAY indicates
%   the full path to the illumination correction measurment file. If
%   STRINGARRAY has only one element, this will be applied to all Z-Planes
%   of IMAGE. If there is the same number of elements in STRINGARRAY as
%   there are Z-planes, they are applied per z-plane (assuming that index
%   of z-plane corresponds to index within STRINGARRAY). Note that the
%   usage of STRINGARRAYs will lead to a reloading of the reference images
%   from the file, slowing down computational time at the benefit of
%   temporarily requiring less RAM (which can be a serious problem when
%   working with large 3D images)
%
%   H = APPLYNBBSILLUMCORRCP3D(Image,MeanImage,StdImage) corrects with the
%   reference MEANIMAGE and reference STDIMAGE each z-plane of IMAGE.
%
%   H = APPLYNBBSILLUMCORRCP3D(Image,cellMeanImage,cellStdImage) corrects
%   each zPlane of Image with another MEANIMAGE and STDIMAGE provided in
%   CELLMEANIMAGE and CELLSTDIMAGE. It assumes that z-plane number equals
%   index within respective cells.
%   ----------------------------------
%   Information about module:
%   Performes NBBS Illumination correction with 2/3D Images. Part of CP3D,
%   adaptation to 3D by TS




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 2
        if iscellstr(varargin)
            bnLoadImages = true;
            CorrPath = varargin;
            numCorrection = length(varargin);
            if min(size(varargin)>1) || ndims(varargin)>2;
                error('Input Filenames must be onedimensional stringarray');
            end
        else
            error('Input must contain either mean and stdev images or Filepath');
        end
    case 3
        if xor(iscell(varargin{1}),iscell(varargin{2}))==true
            error('Mean and Stdev must be provided in same format');
        elseif iscellstr(varargin{1})
            error('Textarray and reference images must not be mixed.');
        elseif iscell(varargin{1}) && iscell(varargin{2})
            if size(varargin{1}) ~= size(varargin{2})
                error('Cells containing reference mean and std are not of equal size.');
            elseif min(size(varargin{1}))>1
                error('Cells containing reference mean and std must be linear.');
            else
                numCorrection = max(size(varargin{1}));
                bnLoadImages = false;
            end
        elseif isnumeric(varargin{1}) && isnumeric(varargin{2})
            if size(varargin{1}) ~= size(varargin{2})
                error('Cells containing reference mean and std are not of equal size.');
            else
                numCorrection = 1;
                bnLoadImages = false;
            end
        else
            error('Could not determine input format of reference for correction.');
        end
    case 1
        error('No reference for correction specified.');
    otherwise
        error('Too many input variables');
end

% Correspond layers in input image with correction references.
numPlanes = size(Image,3);

if numPlanes == numCorrection
    bnCorrectIndividually = true;
elseif numCorrection == 1;
    bnCorrectIndividually = false;
else
    error('There must either be one correction for all z-planes ore one correction per each z-plane. Please check input');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PERFOM ILLUMINATION CORRECTION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preprocess Image
Image = double(Image);
H = Image;      % create copy of input image matrx, where 0 entries will become 1 not that this is important for the log operation later
H(H==0) = 1;

% Load correction for first z-plane
switch bnLoadImages
    case true
        load(CorrPath{1});
        CorrMean = stat_values.mean;
        CorrStd = stat_values.std;
        clear stat_values;
    case false
        if numCorrection == 1;
            CorrMean = varargin{1};
            CorrStd = varargin{2};
        else
            CorrMean=varargin{1}{1};
            CorrStd=varargin{2}{1};
        end
end

% Perform correction
for k=1:numPlanes
    if k~=1 && bnCorrectIndividually == true    % if stacks should be corrected individually, update correction matrices
        switch bnLoadImages     % Determine if Images have to be reloaded from harddrive or are part of input variable
            case true
                load(CorrPath{k});             
                CorrMean = stat_values.mean;
                CorrStd = stat_values.std;
                clear stat_values;
            case false
                CorrMean=varargin{1}{k};
                CorrStd=varargin{2}{k};
        end
    end
        
    H(:,:,k) = (log10(H(:,:,k))-CorrMean)./CorrStd;
    H(:,:,k) = (H(:,:,k).*mean(CorrStd(:)))+mean(CorrMean(:));
end

H = 10.^H;


end