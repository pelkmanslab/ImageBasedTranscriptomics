% Results of individual jobs of Cell Profiler have to be joined together.
% The format of this output will depend upon your job manager. Currently
% only support for iBrain is implemented. 
%
% If a custom spot bias correction is desired, please change
% VisualizePosition.m in such a way that the data can be read from the 
% input variable "XorFile"
%
% Spot positions are encoded by iBrain according to the following scheme:
% Filename: 
% ['Measurments' ObjectName num2str(IncrementOfThreshold) '_Location.mat' ]
% where OBJECTNAME is a string of 'ScanSpot' and an identifier of the plate
% e.g.: ScanSpotCP102r1ab
% and where INCREMENTOFTHRESHOLD is an integer number describing, which
% threshold was tested, eg. 1 for first threshold, 2 for second , ....
%
% This file contains a structure called
% handles
% 
% handles.Measurements.(OBJECTNAME) where 
% handles.Measurements.(OBJECTNAME).Location is a cell where each field contains the output of individual CP cycles 
% as a matrix of objectnumber x 2.
% The first column is the X coordinate of each spot, whereas the second is the Y
% 
% In addition
% handles.Measurements.(OBJECTNAME).LocationFeatures is a cell describing the measurements
% Its content is {'CenterX','CenterY';};
% 
% Please see the example file 
% Measurements_ScanSpotCP102r1ab1_Location as a guideline. Note that
% it only has the spot bias at a single threshold
%
% Useage of iBrain is recommended. For more information about iBRAIN and
% how to obtain it, visit http://www.imls.uzh.ch/research/pelkmans.html