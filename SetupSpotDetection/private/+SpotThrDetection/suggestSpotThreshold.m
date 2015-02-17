function numSpotThreshold = suggestSpotThreshold(ObjectCount,TestedThresholds, MaxSpots, fractionBelowMaxSpots,ddCount,fractionddCount,SliderSize)
% NUMSPOTTHRESHOLD is the value of the lowest Threshold that fullfills the
% criteria for a good threshold. The criteria are 
%   x low amount of spots
%   x stability of amount of spots against variations of threshold



numImageGroups = length(ObjectCount);
numTestedThresholds = length(TestedThresholds);

% Initialize all image groups as permitted and exclude. note that this
% allows to ignore image groups for which no threshold was specified (set
% to NaN).
isAllowed = true(numImageGroups,numTestedThresholds);

% Below MaxSpots Threshold
detMaxSpotseThrGroupIX = find(~isnan(MaxSpots));
if any(detMaxSpotseThrGroupIX)
    for k=1:length(detMaxSpotseThrGroupIX)
        CountIsBelow = ObjectCount{detMaxSpotseThrGroupIX(k)} < MaxSpots(detMaxSpotseThrGroupIX(k));
        fractionBelow = sum(CountIsBelow,1) ./ size(CountIsBelow,1);
        permittedFraction = fractionBelowMaxSpots(detMaxSpotseThrGroupIX(k));
        SetIsPermitted = fractionBelow > permittedFraction;
        % use last threshold below maximum of allowed spots. If thresholds are
        % very low it is possible to have a huge continuous object, filling up
        % the image. in this case the count would be below and before a peak
        % of increasing object numbers.
        numLastAboveThreshold = find(~SetIsPermitted,1,'last');  % instead of checking for first value below upper limit of neg control, determine the last one and add one to indicate the first which after the last non-permitted
        if numLastAboveThreshold == numTestedThresholds  % if no permitted threshold
            isAllowed(detMaxSpotseThrGroupIX(k),:) = false;
        else   % If permissive threshold was found
            isAllowed(detMaxSpotseThrGroupIX(k),1:numLastAboveThreshold) = false;
        end
    end
end

% Change of change of object count (how stable)
detddCountGroupIX = find(~isnan(ddCount));
if any(detddCountGroupIX)
    for k=1:length(detddCountGroupIX)
        % determine change of change of object count and whether fraction
        % of image groups is below limat at a given threshold
        ddObjectCounts = diff(ObjectCount{detddCountGroupIX(k)},2,2);
        plot(abs(ddObjectCounts)');
        figure;
        sdd = sum(abs(ddObjectCounts),2);
        plot(sdd);
        isBelow = ddObjectCounts < ddCount(detddCountGroupIX(k));
        fractionBelow = sum(isBelow,1) ./ size(isBelow,1);
        permittedFraction = fractionddCount(detddCountGroupIX(k));
        SetIsPermitted = fractionBelow > permittedFraction;
        
        % Determine where change is robust over a given slider Size
        numddObjects = size(ddObjectCounts,2);
        isInFrontOfPermittedRange = false(1,numddObjects);
        for r=1:(numddObjects-SliderSize);
            if sum(SetIsPermitted(r:(r+SliderSize-1)))==SliderSize;
                isInFrontOfPermittedRange(1,r)=true;
            end
        end
        
        % update list of permitted thresholds
        isAllowedByddObjects = [false isInFrontOfPermittedRange false]; % note that each difference removes one data point
        isAllowed(detddCountGroupIX(k),:) = isAllowed(detddCountGroupIX(k),:) | isAllowedByddObjects;
    end
end

% find Threshold which fullfils criteria for all image sets
AllowedEveryWhere = sum(isAllowed,1) == numImageGroups;
if any(AllowedEveryWhere)
    numSpotThreshold = TestedThresholds(find(AllowedEveryWhere,1,'first'));
else
    numSpotThreshold = nan;
end
end