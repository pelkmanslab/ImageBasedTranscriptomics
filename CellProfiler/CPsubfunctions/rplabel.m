function imLabel = rplabel(imLogical, imIntensity, Property, logarithm)

%REGIONPROPS_LABEL_IMAGE creates a label image based on a measurement of
%the regionprops function.
%
%   When PROPERTY is provided as string:
%   IMLABEL = REGIONPROPS_LABEL_IMAGE(IMLOGICAL, IMINTENSITY, PROPERTY, LOGARITHM) calls the
%   regionprops function using the input images IMLOGICAL and IMINTENSITY and the input
%   property PROPERTY. It returns a matrix IMLABEL, of the same size as IMLOGICAL,
%   containing labels of the measured PROPERTY for the connected objects in IMLOGICAL. 
%   
%   When PROPERTY is provided as matrix:
%   IMLABEL = REGIONPROPS_LABEL_IMAGE(IMLOGICAL, IMINTENSITY, PROPERTY, LOGARITHM)
%   returns matrix IMLABEL whithout calling the regionprops function.
%   
%   Input: 
%   - imLogical: binary image
%   - imIntensity: intensity image, if you don't want to measure intensities provide
%     empty matrix [] as second input
%   - property: string, e.g. 'Area', 'Eccentricity', 'MeanIntensity', etc.
%     or matrix (when properties were already calculated)
%   - logarithm (optional): string, either 'two' for log2, 'ten' for log10,
%     or 'nat' for log
% 
%   Output:
%   imLabel: label image containing labels of the measured property.
%   (Optionally, output is given in logarithmic form.)


if isempty(imIntensity)
    imIntensity = zeros(size(imLogical));
end

if nargin == 3
    useLog = false;
elseif nargin == 4
    useLog = true;
end
    
if ischar(property)
    matProperty = cell2mat(struct2cell(regionprops(imLogical,imIntensity,property)))';
elseif ismatrix(property)
    matProperty = property;
end
imLabel = bwlabel(imLogical);
Index = unique(imLabel);
Index(Index==0) = [];
for t = 1:length(Index)
    if useLog
        if strcmp(logarithm,'two')
            imLabel(imLabel==Index(t)) = log2(matProperty(t));
        elseif strcmp(logarithm,'ten')
            imLabel(imLabel==Index(t)) = log10(matProperty(t));
        elseif strcmp(logarithm,'nat')
            imLabel(imLabel==Index(t)) = log(matProperty(t));
        end
    else
        imLabel(imLabel==Index(t)) = matProperty(t);
    end
end
