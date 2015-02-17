P{1}.name = 'MyExamplePlate';     % optionally one might have traced a specific plate identifier in the Pipeline setting file and its outcome
P{1}.strBatch = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleSpotDetectionPipeline\BATCH');    % name of the folder which contails the localization of spots
cellObjName = {'UncorrectedSpots';'NotdeblendedSpots';};   % name of the different spot identifications (see PreCluster_Exp_IdentifySpots pipeline)
ImageDimensions = [2160 2560];       % [rows columns] of single images obtained by your microscope (pixel)


%%%%%%%%%%%%%%%%%%%%%%%%
warning('Note that this only is an example to illustrate the code. In order to work in practice, you must learn the correction function from approx. 10,000 different images.')
%%%%%%%%%%%%%%%%%%%%%%%%


outDir = nnpc('\\195.176.109.11\biol_uzh_pelkmans_s7\Data\Users\RNAFish\MethodsPaper\ExampleDataSet\ExampleOutput\BiasCorrection');
ensurePresenceOfDirectory(outDir);

numDifferentSpotClassifications = length(cellObjName);

for j=1:length(P);
        CurrPlate = P{j}.name ;
        cellDescription = {CurrPlate;CurrPlate;CurrPlate;CurrPlate};
        cellBatchPath = repmat({P{j}.strBatch},[2,1]);   % note that it would be possible to unite spots from different BATCH folders by changing this line
        SpotThrDetection.visualizeMultipleSpotBiases(cellBatchPath,cellObjName,cellDescription,ImageDimensions)
        gcf2pdf(outDir,CurrPlate);      %  make pdf
end

%%%%%%%%%%%%%%%%%
warning('Note that this only is an example to illustrate the code. In order to work in practice, you must learn the correction function from approx. 10,000 different images.')
%%%%%%%%%%%%%%%%%%%%
