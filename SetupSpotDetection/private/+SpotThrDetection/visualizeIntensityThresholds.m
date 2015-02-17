function visualizeIntensityThresholds(minIntensity, maxIntensity, vIntensityBoundaries,XLimits)
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html

if nargin < 4
    bnCustomXLimits = false;
else
    bnCustomXLimits = true;
    
end

if size(minIntensity,2) ~= size(maxIntensity,2)
    error('The number of image groups within minIntensity and maxIntensiy is different');
end

numGraphs = size(minIntensity,2)*2;

% get lowest and highest intensity values
absMin = min(...
    [cell2mat(cellfun(@(x) min(x), maxIntensity, 'UniformOutput', false))...
    cell2mat(cellfun(@(x) min(x), minIntensity, 'UniformOutput', false))]);
absMax = max(...
    [cell2mat(cellfun(@(x) max(x), maxIntensity, 'UniformOutput', false))...
    cell2mat(cellfun(@(x) max(x), minIntensity, 'UniformOutput', false))]);
tmpExtend = floor((absMax-absMin)*0.05);
absMin = absMin-tmpExtend;
absMax = absMax+tmpExtend;


figure;

for k=1:size(minIntensity,2)
    subplot(numGraphs,1,k*2-1)
    [n,xout]=hist(minIntensity{k},absMin:2:absMax);
    bar(xout,n);
    maxn = max(n);
    line([vIntensityBoundaries(1),vIntensityBoundaries(1)],[0,maxn],'Color','g','Linewidth',3);
    line([vIntensityBoundaries(2),vIntensityBoundaries(2)],[0,maxn],'Color','r','Linewidth',3);
    line([vIntensityBoundaries(3),vIntensityBoundaries(3)],[0,maxn],'Color','g','Linewidth',3);
    line([vIntensityBoundaries(4),vIntensityBoundaries(4)],[0,maxn],'Color','r','Linewidth',3);
    xlabel(['Image Minima of Image set ' num2str(k)]);
    
    if bnCustomXLimits == true
        xlim(XLimits);
    end
    
    if k==1
        strIntensityBoundaries = num2str(vIntensityBoundaries);
        title(['IntensityThresholds: ' strIntensityBoundaries]);
        
    end
    
end

for k=1:size(maxIntensity,2)
    subplot(numGraphs,1,k*2)
    [n,xout]= hist(maxIntensity{k},absMin:2:absMax);
    bar(xout,n);
    maxn = max(n);
    line([vIntensityBoundaries(1),vIntensityBoundaries(1)],[0,maxn],'Color','g','Linewidth',3);
    line([vIntensityBoundaries(2),vIntensityBoundaries(2)],[0,maxn],'Color','r','Linewidth',3);
    line([vIntensityBoundaries(3),vIntensityBoundaries(3)],[0,maxn],'Color','g','Linewidth',3);
    line([vIntensityBoundaries(4),vIntensityBoundaries(4)],[0,maxn],'Color','r','Linewidth',3);
    xlabel(['Image Maxima of Image set ' num2str(k)]);
    
    if bnCustomXLimits == true
        xlim(XLimits);
    end
    
end

end