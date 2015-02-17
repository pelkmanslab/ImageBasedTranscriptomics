function [RescaleThr] = suggestImageIntensityBoundaries(A,qA,B,qB,C,qC,D,qD)
% creates RESCALETHR for OBJBYFILTER(). Based upon the quantiles QX for
% each X. The reason for quantiles is to prevent highly unusual images to affect
% global intensity Threhsolds.
%
% individual X can also be NaN. Then RESCALETHR will have NaN at this
% position, whereby this specific Threshold, but not the non-nan elements,
% will be ignored by OBJBYFILTER

RescaleThr = nan(1,4); % initialize

if ~(isempty(qA) || isnan(qA))
  RescaleThr(1) = quantile(A,qA);  
end

if ~(isempty(qB) || isnan(qB))
  RescaleThr(2) = quantile(B,qB);  
end

if ~(isempty(qC) || isnan(qC))
  RescaleThr(3) = quantile(C,qC);  
end

if ~(isempty(qD) || isnan(qD))
  RescaleThr(4) = quantile(D,qD);  
end

end
