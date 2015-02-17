function [intRow, intColumn, intImagePosition, intTimepoint, intZstackNumber, intChannelNumber, strMicroscopeType, strWellName, intActionNumber] = MetaFromImageName(strImageName)
% METAFROMIMAGES collects metainformation about image acquistion from the
% file name. It serves as a hub for multiple different functions,
% previously developed within the lab, where features are derived from
% parsing the file names. STRIMAGENAME is the name of the file. You might
% add custom modifications to extract your metainformation of interest.
%
%%%%%%%% IMPORTANT: MAKE YOUR CUSTOM ADJUSTMENTS %%%%%%%%%%%%%%
% a) Use regular expressions (see matlab help) to obtain the information. 
% This might look something like
% strChannelMatch = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})L(\d{2})A(\d{2})Z(\d{2})C(\d{2})', 'Tokens');
% b) Some of the output arguments are actually not required for Spot detection, 
% so you can set them to default value (indicated by * in
% description)
%
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html
% *************************************************************************
%
%   NAME                TYPE       
%   intChannelNumber    double      
%   Number describing the color, eg. 1 for blue, 2 for green, 3 for red, 4 
%   for far red,.... (might be custom)
%
%   intZstackNumber     double      
%   Number describing the Z-plane. NaN if not suitable, otherwise: start 
%   with 1 and increment by 1 corresponding to subsequent stage positions.
%   *(1) *only if no 3D analysis
%   
%   intActionNumber     double      
%   Number describing the action. eg.: if multiple channels were acquired 
%   at the same time by parallel optics/cameras. NaN if not suitable. *(1)
%
%   intImagePosition    double
%   Number describing the acquisition site within a single well. 
%
%   strMicroscopeType   character
%   Name of microscope. Mostly used for better error messages. *('myMic')
%
%   intRow      double
%   Row of well. Row A is 1. Row B is 2. Row C is 3.....
%
%   intColumn       double
%   Column of well.
%
%   strWellName   double
%   Name of Well accoding with combination of letter (for rows) and number
%   (for columns) eg. A05 . *('A01')
%
%   int Timepoint   double *
%   Number indicating the timepoint of a time series. Starts with 1 and
%   increments by 1 for subsequent time points. *(1)


% Add security and informative error message for publication
if isempty(which('check_image_channel')) || isempty(which('check_image_channel')) || isempty(which('filterimagenamedata'))
    error('Sorry. You have to do some coding: You have to write custom functions, which obtain the metadata (such as image channel, Position in Z-stack...) from the filename (or alternative source). Place them in MetaFromImageName.m')
else
    % Call functions, which anlyse image file names for metadata
    [intChannelNumber,intZstackNumber,intActionNumber] = check_image_channel(strImageName);
    [intImagePosition,strMicroscopeType] = check_image_position(strImageName);
    [intRow, intColumn, strWellName, intTimepoint] = filterimagenamedata(strImageName);
end

end