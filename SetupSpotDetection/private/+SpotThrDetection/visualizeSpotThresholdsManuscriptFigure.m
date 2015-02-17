function visualizeSpotThresholdsManuscriptFigure(ObjectCount, TestedThresholds, strMethod, numChosenThreshold,ZoomOnY,customxLimits)

numGraphs = 1 + length(ZoomOnY);
numGroups = length(ObjectCount);
xLimits = [TestedThresholds(1) TestedThresholds(end)];
if ~isempty(customxLimits)
    xLimits = customxLimits;
end

% Colourmap, note that first and second row are blue and red respectively;
% others correspond to color palette jet
Colours = [0.5,0.5,0.5;1,0,0;0,0,0.687500000000000;0,0,0.750000000000000;0,0,0.812500000000000;0,0,0.875000000000000;0,0,0.937500000000000;0,0,1;0,0.0625000000000000,1;0,0.125000000000000,1;0,0.187500000000000,1;0,0.250000000000000,1;0,0.312500000000000,1;0,0.375000000000000,1;0,0.437500000000000,1;0,0.500000000000000,1;0,0.562500000000000,1;0,0.625000000000000,1;0,0.687500000000000,1;0,0.750000000000000,1;0,0.812500000000000,1;0,0.875000000000000,1;0,0.937500000000000,1;0,1,1;0.0625000000000000,1,0.937500000000000;0.125000000000000,1,0.875000000000000;0.187500000000000,1,0.812500000000000;0.250000000000000,1,0.750000000000000;0.312500000000000,1,0.687500000000000;0.375000000000000,1,0.625000000000000;0.437500000000000,1,0.562500000000000;0.500000000000000,1,0.500000000000000;0.562500000000000,1,0.437500000000000;0.625000000000000,1,0.375000000000000;0.687500000000000,1,0.312500000000000;0.750000000000000,1,0.250000000000000;0.812500000000000,1,0.187500000000000;0.875000000000000,1,0.125000000000000;0.937500000000000,1,0.0625000000000000;1,1,0;1,0.937500000000000,0;1,0.875000000000000,0;1,0.812500000000000,0;1,0.750000000000000,0;1,0.687500000000000,0;1,0.625000000000000,0;1,0.562500000000000,0;1,0.500000000000000,0;1,0.437500000000000,0;1,0.375000000000000,0;1,0.312500000000000,0;1,0.250000000000000,0;1,0.187500000000000,0;1,0.125000000000000,0;1,0.0625000000000000,0;1,0,0;0.937500000000000,0,0;0.875000000000000,0,0;0.812500000000000,0,0;0.750000000000000,0,0;0.687500000000000,0,0;0.625000000000000,0,0;0.562500000000000,0,0;0.500000000000000,0,0;];

% expand colour palette if necessary
lColours = size(Colours,1);
if numGroups > lColours;
    j = ceil(numGroups/lColours);
    Colours = repmat(Colours,j,1);
end

% see if permitted Threshold is present
if ~isnan(numChosenThreshold) || ~isempty(numChosenThreshold);
    bnThresholdIsOk = true;
    if numChosenThreshold > xLimits(2)
        xLimits(2) = numChosenThreshold;
    end
else
    bnThresholdIsOk = false;
end


% Visualize
switch strMethod
    % Display absolute count of Spots
    case 'Absolute'
        figure;
        for k=1:numGraphs
            maxValue = 0;
            for l=1:numGroups
                subplot(numGraphs,1,k);
                plot(TestedThresholds,ObjectCount{l},'Color',Colours(l,:));
                hold on;
                if k==1;
                    tmpMax = max(ObjectCount{l}(:));
                    maxValue = max(maxValue, tmpMax);
                    title(['Number of Spots among individual images vs. threshold of detection: ' num2str(numChosenThreshold)]);
                end
                
                     MedianCount = median(ObjectCount{l},1);
                    plot(TestedThresholds,MedianCount,'b','Linewidth',3);

                
            end
            
            
            if k>1  % while the first subgraph shows the complete range, the others have a user defined one
                maxValue = ZoomOnY(k-1);
            end
            
            % display chosen threshold and adapt x-scale, if neceessary
            if bnThresholdIsOk == true
                line([numChosenThreshold,numChosenThreshold],[0 maxValue],'Color','g','Linewidth',3); % indicate recommeneded threshold
            end
            
            xlim(xLimits); % force all graphs to show full range of tested parameters
            ylim([0 maxValue]);
        end
        
    case 'Relative'
        if bnThresholdIsOk == false
            fprintf('Visualization of Spot count relative to threshold is not possible because no threshold has been identified \n');
        else
            % Find Threshold closest to chosen one
            [~, PositionOfThreshold] = min(abs(TestedThresholds- numChosenThreshold));
            
            figure;
            for k=1:numGraphs
                maxValue = 0;
                for l=1:numGroups
                    tmp = repmat(ObjectCount{l}(:,PositionOfThreshold),1,size(TestedThresholds,2));
                    RelativeCount = ObjectCount{l}./tmp;
                    subplot(numGraphs,1,k);
                    plot(TestedThresholds,RelativeCount,'Color',Colours(l,:));
                    hold on;
                    MedianCount = median(RelativeCount,1);
                    plot(TestedThresholds,MedianCount,'b','Linewidth',3);

                    
                    if k==1;
                        tmpMax = max(RelativeCount(:));
                        maxValue = max(maxValue, tmpMax);
                        title('Number of Spots among individual images vs. threshold of detection.');
                    end
                end
                
                if k>1  % while the first subgraph shows the complete range, the others have a user defined one
                    maxValue = ZoomOnY(k-1);
                end
                
                line([numChosenThreshold,numChosenThreshold],[0 maxValue],'Color','k','Linewidth',3);
                
                xlim(xLimits); % force all graphs to show full range of tested parameters
                ylim([0 maxValue]);
            end
        end
end

end