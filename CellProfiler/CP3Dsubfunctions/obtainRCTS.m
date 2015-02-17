function [subRCTS, subMax, intZPos] = obtainRCTS(FileList)
% Support function for CP3D
% For image names provided in the input array File list, this function
% creates a linear index consisting of RowColumnTimepointSite, which can be 
% converted back by subMax. The third output indicates the z plane.
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html
% *************************************************************************

[intRow, intColumn, intImagePosition, intTimepoint, intZstackNumber, ~, ~, ~, ~] =  cellfun(@MetaFromImageName,FileList,'UniformOutput',false);
intRefRCTS = [cell2mat(intRow) cell2mat(intColumn) cell2mat(intTimepoint) cell2mat(intImagePosition)];
subMax = max(intRefRCTS,[],1);
subRCTS = sub2ind2(subMax,intRefRCTS);
intZPos = cell2mat(intZstackNumber);
end