function handles = MeasureSingerFeatures(handles)

% Help for the MeasureSingerFeatures
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures features described in Park et al. 2012 Cell Reports, which deals
% with RNA Localization
% *************************************************************************
%
% see Park et al. 2012 Cell Reports
%
% Measurement:             Feature Number:
% PolarizationIndex                 |   1
% SecondMomentOfChild               |   2
% SecondMomentOfHypotheticalUniform    |   3
% DispersionIndex           |   4
%
% POLARIZATION INDEX
% Measure of the distance of the centroid of all children compared to the
% centroid of the cell. This is normalized by the root-mean-square distance
% of all pixels of the parent. Notably this might be misleading in some cases since
% directional information is not used. Thus it could be possbile that
% children, where the difference to centroid might have the same amplitude
% as any children, but the direction of the extension is opposite seem
% unpolarized whereas they might be opposingly polarized.
%
% DISPERSION INDEX
% Measure how evenly distributed individual children are compared to the
% centroid of all children
%
% SECOND MOMENT OF INTENSITIES
% Measure of the distribution of intensities within one object compared to
% the centroid of the object
%
%
% How it works:
% See initial publication for formulas.
%
%
% Authors:
%   TS
%   NB   
%
% $Revision: 4526 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What are the parent objects?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
iParentName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the child objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
iChildName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which image should be used for the intensity based features?
%infotypeVAR03 = imagegroup
iImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu



%%%VariableRevisionNumber = 2


% Set up the window for displaying the results
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    CPfigure(handles,'Text',ThisModuleFigureNumber);
    columns = 1;
end

% Get current Cycle ID
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & DATA HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Retreive Data   %%%%%%%%%%%%%%%%%
%chech if child have parents and find the index.
% [NB] What is the object has two parents! Need to look for the correct parent. 
IXParentObject = find(cell2mat(cellfun(@(x) strcmp(x, iParentName), handles.Measurements.(iChildName).ParentFeatures, 'uniformoutput',false)),1,'first');
if isempty(IXParentObject)
   error('%s: Please relate the clidren objects to the parent objects!!!. Amature!.',mfilename)
end

% Parent ID of individual child objects
ParentOfChild = handles.Measurements.(iChildName).Parent{SetBeingAnalyzed}(:,IXParentObject);

% Centroid of Objects
CentroidOfParent = handles.Measurements.(iParentName).Location{SetBeingAnalyzed};
CentroidOfAllChildren  = handles.Measurements.(iChildName).Location{SetBeingAnalyzed};

% [TS 120906] commented out check for dimensions to allow empty
% segmentation images, also this module is currently only used in 2D
% Pipelines
% if size(CentroidOfParent,2)~=2     % force the use of 2D images
%     error(['Image processing was canceled in the ', ModuleName, ' module. There was a problem with the dimensions. The Parent Object is not 2D.'])
% end
% 
% if size(CentroidOfAllChildren,2)~=2     % force the use of 2D images
%     error(['Image processing was canceled in the ', ModuleName, ' module. There was a problem with the dimensions. The Child Object is not 2D.'])
% end

% Number (Amount) of Parent objects
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,iParentName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {iParentName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
NumParents = handles.Measurements.Image.ObjectCount{SetBeingAnalyzed}(1,column);

%%%%%%%%%%%%%%%% Retreive Images / Segmentation   %%%%%%%%%%%%%%%%%

% Opens the image used for determining intensity features.
OrigImage = CPretrieveimage(handles,iImageName,ModuleName);
ImageDim = size(OrigImage);

% Obtain Segmentation of Parent
LabelMatrixParent = CPretrieveimage(handles,['Segmented', iParentName],ModuleName,'MustBeGray','DontCheckScale');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get properties of Parent Object
props = regionprops(LabelMatrixParent,'PixelList');

% Determine, which parents have children and store this information as
% logical where position corresponds to individual parents
ParentsWithChildren = unique(ParentOfChild);
f = ParentsWithChildren > 0;
ParentsWithChildren = ParentsWithChildren(f); % Backround not accepeted as parent
bnHasChild = false(NumParents,1);
bnHasChild(ParentsWithChildren) = true; clear ParentsWithChildren;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EXPLAINATION OF FURTHER CODE
%%%%% Children are presorted to enable simple looking up without additional
%%%%% loops or redundant search for cognate children; Loop through each
%%%%% parent and obtain required information of this parent and its
%%%%% children, calculate features according to description of Park et al.


%%% Initialize Output

% Initialize Output. If output is n/a, eg. if no child use NaN as output
% since 0 can also be the result of some measurements, such as for a
% hypothetically not polarized child
LabelOfFeatures    = {'PolarizationIndex',...
    'SecondMomentOfChild','SecondMomentOfHypotheticalUniformDist','DispersionIndex','SecondMomentOfIntensities'};
NumOutputs = length(LabelOfFeatures);
SingerMeasurements=NaN(NumParents,NumOutputs);

if NumParents> 0 % If there is at least one parent object
    %%% Preprocess Child data and sort them according to parentid
    % Filter Child data for children with non-background parent
    f = ParentOfChild > 0;
    ParentOfChild = ParentOfChild(f);
    CentroidOfAllChildren = CentroidOfAllChildren(f,:);
    % Sort Child Centroids according to their Parents
    [ParentOfChild, sortIX] = sort(ParentOfChild);
    CentroidOfAllChildren = CentroidOfAllChildren(sortIX,:);  clear sortIX;
    % Obtain positions where there are new parents, where there is last child
    % of a parent
    
    % below line added by berend
    if ~isempty(ParentOfChild)

        [~, IXFirstChild, ~]= unique(ParentOfChild(:,1),'first');
        IXLastChild = IXFirstChild(2:end)-1;
        IXLastChild = [IXLastChild; size(CentroidOfAllChildren,1)];
        j=1; %
        for k=1:NumParents
            %%% SUMMARY
            % PART 1 yields features independent of the presence of a child
            % PART 2 yields features dependent of the presence of a child

            %%% PART 1: INDEPENDENT UPON CHILD

            % PixelList of current object, here in format of [x y], which will be
            % directly applicable for distance measurements
            CoordinatesAllPixelsOfAParent = props(k).PixelList;
            % Convert PixelsIDs into linear index. note that x is only the second
            % dimension for adressing matrices in matlab, but the first dimension
            % of PixelList, therefore the following line swapps
            subCoordinatesAllPixelsOfAParent = sub2ind(ImageDim,CoordinatesAllPixelsOfAParent(:,2),CoordinatesAllPixelsOfAParent(:,1));
            % Centroid of Current Parent
            CentroidCurrentParent = CentroidOfParent(k,:);

            %%% Feature: Second momentum of intensity

            % Obtain Intensities of individual pixels.
            IntensitiesOfAllParentPx = OrigImage(subCoordinatesAllPixelsOfAParent);
            % Total Intensity within the parent object.
            TotalIntensityOfParent = sum(IntensitiesOfAllParentPx);

            % Square distances of individual pixels to centroid of Cell
            tmp = CoordinatesAllPixelsOfAParent - repmat(CentroidCurrentParent,size(CoordinatesAllPixelsOfAParent,1),1);
            tmp = tmp.^2;
            SquareOfDistances = sum(tmp,2);  clear tmp;

            % Second moment of intensities
            tmp = SquareOfDistances .* IntensitiesOfAllParentPx;
            SecondMomentOfIntensites = sum(tmp) ./TotalIntensityOfParent;   clear tmp;

            % Output
            SingerMeasurements(k,5) = SecondMomentOfIntensites;


            %%% PART 2: DEPENDENT UPON CHILD
            if bnHasChild(k) == true
                CentroidOffspring = CentroidOfAllChildren(IXFirstChild(j):IXLastChild(j),:);
                j = j+1;  % increase index within lookup-table

                %%% Feature: Polarization Index

                % Size Of Polarization Vector
                MeanPosChild = mean(CentroidOffspring,1);
                PolarizationVector = (MeanPosChild-CentroidCurrentParent);
                SizeOfPolarizationVector = sqrt(sum(PolarizationVector.^2));

                % Radius of gyration: " The radius of gyration Rgcell is calculated
                % by the root-mean-square distance of all pixels within the cell
                % from the centroid of the cell" (Park et al.)
                tmp = CoordinatesAllPixelsOfAParent - repmat(CentroidCurrentParent,size(CoordinatesAllPixelsOfAParent,1),1);  
                tmp = sum(sum(tmp.^2,2)) * 1/size(CoordinatesAllPixelsOfAParent,1); 
                RadiusOfGyration = sqrt(tmp); clear tmp;

                % PolarizationIndex
                PolarizationIndex = SizeOfPolarizationVector ./ RadiusOfGyration;

                % Output
                SingerMeasurements(k,1) = PolarizationIndex;

                %%% Feature: Dispersion Index

                % Second moment of Child
                tmp = std(CentroidOffspring,1,1);
                SecondMomentOfChild = sum(tmp.^2); clear tmp;

                % Output
                SingerMeasurements(k,2) = SecondMomentOfChild;

                % Second moment of a hypothetical uniform distribution
                tmp = CoordinatesAllPixelsOfAParent - repmat(MeanPosChild,size(CoordinatesAllPixelsOfAParent,1),1);
                SecondMomentOfHypothetical = sum(sum(tmp.^2)) * 1/size(CoordinatesAllPixelsOfAParent,1); clear tmp;      

                % Output
                SingerMeasurements(k,3) = SecondMomentOfHypothetical;

                % Dispersion Index
                DispersionIndex = SecondMomentOfChild./SecondMomentOfHypothetical;

                % Output
                SingerMeasurements(k,4) = DispersionIndex;
            end    
        end
    end
else
    SingerMeasurements = NaN(1,NumOutputs);
end

    
    
    %%% Save measurements
    handles.Measurements.(iParentName).(['Singer_',iImageName,'Features']) = LabelOfFeatures;
    handles.Measurements.(iParentName).(['Singer_',iImageName])(SetBeingAnalyzed) = {SingerMeasurements};
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% REPORT MEASUREMENTS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% This code is a modified copy of the report of the MeasureObjectIntensity module
    
    
    %%% Report measurements
    if any(findobj == ThisModuleFigureNumber);
        FontSize = handles.Preferences.FontSize;
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            delete(findobj('parent',ThisModuleFigureNumber,'string','R'));
            delete(findobj('parent',ThisModuleFigureNumber,'string','G'));
            delete(findobj('parent',ThisModuleFigureNumber,'string','B'));
        end
        %%%% This first block writes the same text several times
        %%% Header
        
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0 0.95 1 0.04],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string',sprintf(['Average intensity features for ', iImageName,', cycle #%d'],handles.Current.SetBeingAnalyzed));
        
        %%% Number of objects
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.85 0.3 0.03],...
            'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string','Number of objects:');
        
        %%% Text for Basic features
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.8 0.3 0.03],...
            'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string','Intensity feature:');
        for k = 1:NumOutputs
            uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.8-0.04*k 0.3 0.03],...
                'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
                'fontsize',FontSize,'string',LabelOfFeatures{k});
        end
        
        %%% The name of the object image
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.9 0.1 0.03],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string',iParentName);
        
        %%% Number of objects
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.85 0.1 0.03],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'string',num2str(NumParents));
        
        if NumParents > 0
            %%% Singer features
            for k = 1:NumOutputs
                uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.8-0.04*k 0.1 0.03],...
                    'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
                    'fontsize',FontSize,'string',sprintf('%0.2f',nanmean(SingerMeasurements(:,k))));
            end
        end
        %%% This variable is used to write results in the correct column
        %%% and to determine the correct window size
        columns = columns + 1;  % This is seems funny since unused, however the same logic is in CP original MeasureObjectIntensity module. might be important for wrapper?
        
        % [NB] plot scater
        subplot(2,2,2);
        plot(SingerMeasurements(:,1),SingerMeasurements(:,4),'k.')
        xlabel('PolarizationIndex')
        ylabel('DispersionIndex')
        
        % [NB] plot hist
        subplot(2,3,4);
        hist(SingerMeasurements(:,1))
        title('PolarizationIndex')
        
        subplot(2,3,5);
        hist(SingerMeasurements(:,4))
        title('DispersionIndex')
        
        subplot(2,3,6);
        hist(SingerMeasurements(:,5))
        title('SecondMomentOfIntensities')
        
       
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % % %%%%%%%%%%%%%%%%%% ARTICLE  %%%%%%%%%%% Park et al. 2012 Cell Reports
% % %     
% % %     An UnbiasedAnalysisMethod to QuantifymRNA Localization Reveals Its Correlation with Cell Motility
% % % 
% % %     Hye Yoon Park1, 2,
% % %     Tatjana Trcek1,
% % %     Amber L. Wells1,
% % %     Jeffrey A. Chao1,
% % %     Robert H. Singer1, 2, Corresponding author contact information, E-mail the corresponding author
% % % 
% % %     1 Department of Anatomy and Structural Biology, Albert Einstein College of Medicine, Bronx, NY 10461, USA
% % %     2 Gruss Lipper Biophotonics Center, Albert Einstein College of Medicine, Bronx, NY 10461, USA
% % % 
% % %     Received 30 August 2011. Revised 22 November 2011. Accepted 23 December 2011. Available online 16 February 2012. Published online: February 16, 2012. 
% % % 
% % %     http://dx.doi.org/10.1016/j.celrep.2011.12.009, How to Cite or Link Using DOI
% % %     Cited by in Scopus (0)
% % % 
% % %     Permissions & Reprints
% % % 
% % % Summary
% % % 
% % % Localization of mRNA is a critical mechanism used by a large fraction of transcripts to restrict its translation to specific cellular regions. Although current high-resolution imaging techniques provide ample information, the analysismethods for localization have either been qualitative or employed quantification in nonrandomly selected regions of interest. Here, we describe an analytical method for objective quantification of mRNA localization using a combination of two characteristics of its molecular distribution, polarization and dispersion. The validity of the method is demonstrated using single-molecule FISH images of budding yeast and fibroblasts. Live-cell analysis of endogenous ?-actin mRNA in mouse fibroblasts reveals that mRNA polarization has a half-life of ?16 min and is cross-correlated with directed cell migration. This novel approach provides insights into the dynamic regulation of mRNA localization and its physiological roles.
% % % Graphical Abstract
% % % 
% % % Full-size image (29K)
% % % 
% % % Highlights
% % % 
% % % 
% % % Cells achieve polarity in part by localization of mRNA, which allows protein synthesis to be confined to a subcellular compartment (Meignin and Davis, 2010). In the studies of mRNA localization in somatic cells, visualization has played a crucial role since the first observation almost 30 years ago ( [Jeffery et al., 1983] and [Lawrence and Singer, 1986]). Although in situ hybridization is still widely considered as the standard tool, a variety of techniques in imaging and labeling have enabled detection of single RNA molecules not only in fixed cells (Femino et al., 1998) but also in live cells in real time (Bertrand et al., 1998).
% % % 
% % % Although imaging techniques are highly sophisticated, analysis of mRNA localization has been mostly limited by the qualitative interpretation. For instance, transcripts expressed during Drosophila embryogenesis are classified into ?35 localization categories (Lécuyer et al., 2007). A conventional method to quantify RNA localization is based on manual counting of the cells by two independent observers who are blind to the experimental conditions. Although it may eliminate potential bias in data selection and processing, ambiguities still remain due to individual variation. In addition to the need for an objective analysis, the importance of quantitative measurement has been recognized for several reasons. First, localization can be described as a continuous process rather than an “all-or-nothing” occurrence (Zimyanin et al., 2008). Second, a metric for localization could identify distinct populations that may be missed by a binary analysis. Moreover, quantification could be automated to facilitate high-throughput image analysis where differential response of cells could be examined under manifold conditions.
% % % 
% % % A few quantification methods have been suggested in the literature to analyze subcellular localization of mRNA. In one approach, the most dense and least dense regions of the cell were selected, and the ratio of the highest RNA density to the lowest density was calculated (Lawrence and Singer, 1986). Latham et al. (1994) counted a cell as localizing if 80% signal was concentrated in leading lamella area comprising a quarter of the cell area. More recently, Yamagishi et al. (2009) selected regions of interest (ROI), and calculated the ratio of the mRNA concentration at the leading edge and in the perinuclear region. For the transcripts that localize in the vicinity of a certain cellular compartment, the distance between the mRNA and the cellular objects may be used to quantify localization (Jourdren et al., 2010). These varying analysismethods could lead to significantly different conclusions, and demonstrate the need for an objective quantitative analysis of RNA localization.
% % % 
% % % Here we demonstrate an analytical method applied to three different cell types: budding yeast, chicken embryonic fibroblasts (CEF), and mouse embryonic fibroblasts (MEF). We introduce two measures to characterize mRNA distribution, namely polarization and dispersion. When mRNA distribution is asymmetric in a cell, the centroid of the mRNA should deviate from the centroid of the cell. We define the displacement vector pointing from the center of the cell to the center of mRNA as the RNA polarization vector (Figure 1A). The polarization index (PI) is determined by dividing the size of the polarization vector by the radius of gyration of the cell (Rgcell)
% % % View the MathML source
% % % where View the MathML source, View the MathML source, View the MathML source are the coordinates of the centroid for the RNA, and View the MathML source, View the MathML source, View the MathML source are those for the cell in a general three-dimensional case. The radius of gyration Rgcell is calculated by the root-mean-square distance of all pixels within the cell from the centroid of the cell. The displacement between the two centroids is divided by Rgcell in order to assess the polarization normalized to the size and the elongation of the cell. The second quantity, the dispersion of mRNA, is measured by calculating the second moment ?2 of RNA positions:
% % % View the MathML source
% % % where N is the total number of mRNA molecules, xi, yi, zi are the coordinates of the ith mRNA, and View the MathML source, View the MathML source, View the MathML source are the variances of the positions. The second moment is dependent on the shape and size of the cell as well as the RNA distribution. To normalize the effect from the cell morphology, we divide the second moment of mRNA?2 by the second moment of the hypothetical uniform distribution View the MathML source. A binary mask image of each cell is generated and View the MathML sourceis calculated as the second moment of each pixel's coordinates within the mask
% % % View the MathML source
% % % where M is the total number of pixels within the mask, and Xi, Yi, Zi are the coordinates of the ith pixel. Then we define the dispersion index (DI) as
% % % View the MathML source
% % % 
% % % Full-size image
% % % 
% % %     Figure 1. Quantification of ASH1mRNA Distributions in Wild-Type and ?SHE2 Budding Yeasts(A) A schematic showing the polarization vector of mRNA distribution pointing from the centroid of the cell to the centroid of the single mRNA positions.(B) Overlay image of ASH1mRNA molecules (magenta), and nuclei (blue) in a wild-type cell. ASH1 mRNAs localize to the daughter bud tip. Scale bar represents 1 ?m.(C) Scatter plot of polarization index in x axis and dispersion index in y axis for wild-type cells.(D) Overlay image of ASH1mRNA molecules, and nuclei in a ?SHE2 cell. Deletion of She2p causes a complete delocalization of ASH1 transcripts. Scale bar represents 1 ?m.(E) Scatter plot of polarization index versus dispersion index for ?SHE2 cells. The red arrow heads indicate the data points for the cells shown in (B) and (D).See also Figure S1.
% % % 
% % % DI has by definition a value of 1 if the cell has an absolutely uniform distribution of mRNA (Figure S1A). When mRNA is concentrated in a certain region, DI is less than 1 (Figures S1B–S1D). If mRNA is spread around the rim of the cell, DI becomes larger than 1 (Figure S1E).
% % % 
% % % Full-size image
% % % 
% % %     Figure S1. Monte Carlo Simulation of Localization Patterns, Related to Figure 1(A-E) Simulated localization patterns of 1,000 particles in a bounded circle. (F) Each data point indicates the polarization index and the dispersion index of the corresponding localization pattern shown in (A)-(E).
% % % 
% % % To test the two-quantity approach (PI and DI), we analyzed localization of the budding yeast ASH1mRNA. ASH1 expression is cell-cycle regulated and reaches a peak during mitosis. In wild-type (WT) cells, the majority of ASH1 transcripts are actively localized to the bud (Long et al., 1997) (Figure 1B). Among several proteins that regulate ASH1 localization, She2p is the primary RNA-binding protein (Niessing et al., 2004). In the absence of She2p, ASH1mRNA becomes homogeneously distributed between the mother cell and the bud tip (Long et al., 1997) (Figure 1D). We performed single-molecule fluorescence in situ hybridization (FISH) (Femino et al., 1998) on WT and ?SHE2 cells as described previously ( [Trcek et al., 2012] and [Zenklusen et al., 2008]). Multiple fields were imaged and maximum intensity projections of Z stacks were generated for two-dimensional analysis. By using a least-squares Gaussian fitting routine, we obtained the intensity and spatial information of each fluorescent particle and processed them to quantify polarization and dispersion of RNA distribution. In order to analyze the posttranscriptional localization of mRNA, the multiple copies of mRNA at the transcription sites were excluded from the analysis. In wild-type strain, PI was 0.71 ± 0.04 (SEM) and DI was 0.50 ± 0.04. In ?SHE2 strain, PI was 0.16 ± 0.02 and DI was 0.82 ± 0.02. Consistent with the human perception of the representative images (Figures 1B and 1D), the objective metrics indicate that the highly polarized localization of ASH1mRNA in wild-type cells was disrupted by the deletion of SHE2 gene. In order to assess the relationship between the polarization and dispersion of mRNA in each strain, we computed the correlation coefficient of the two metrics. The Pearson correlation coefficient was ?0.72 in wild-type and ?0.13 in the ?SHE2 strain. In wild-type cells, higher polarization of ASH1mRNA was strongly associated with tighter confinement of the mRNA. This negative correlation between the polarization and the dispersion disappeared in ?SHE2 strain (Figures 1C and 1E). The correlation coefficient describes the characteristics of RNA distribution in a certain population of cells.
% % % 
% % % We next applied the analysismethod for the distribution of mRNA in chicken embryonic fibroblasts. By comparing the distribution of ?-actin mRNA and GAPDH mRNA in the same cell, we investigated if localization is a general effect for any mRNA or a specific effect for ?-actin mRNA. Figure 2A shows a representative image of primary CEF cells in which ?-actin mRNA is localized to the leading edge (indicated with arrowheads), whereas GAPDH mRNA is more uniformly distributed. Because the densities of these mRNA species are too high to distinguish individual molecules, the distribution of pixel intensity values were analyzed instead of the distribution of discrete mRNA particles. After background subtraction, the second moment ?2 was calculated by
% % % View the MathML source
% % % where rij is the distance from the centroid of the cell to the pixel (i,j) within the cell boundary, and Iij is the intensity value of the pixel in a two-dimensional image. A representative cell shown in the left side of Figure 2A exhibits polarized localization of ?-actin mRNA (PI = 1.14, DI = 0.35) and more uniform distribution of GAPDH mRNA (PI = 0.33, DI = 1.17). In a randomly-chosen population of CEF cells (n = 99), the mean PI was 0.47 ± 0.03 for ?-actin mRNA and 0.29 ± 0.02 for GAPDH mRNA. In 82% of the cells (in the upper triangle above the red dashed line in Figure 2B), the distribution of ?-actin mRNA was more polarized than GAPDH mRNA. However, the difference in the mean DI was not as significant: 0.98 ± 0.05 for ?-actin mRNA and 0.82 ± 0.04 for GAPDH mRNA. Unlike yeast ASH1mRNA, the mean values of PI and DI for ?-actin mRNA indicate a moderate polarization and a loose dispersion on average because of the intrinsic cell-to-cell variation in primary CEF culture. Previously, Latham et al. (1994) reported that only 30%–35% of primary CEF localized ?-actin mRNA to the cell periphery using a binary counting method. These cells showing ?-actin mRNA localization may represent a subpopulation of motile cells (Kislauskis et al., 1997). Although the cell population was heterogeneous, the correlation coefficient between PI and DI for ?-actin mRNA was highly negative (r = ?0.61) in contrast to the low value for GAPDH mRNA (r = ?0.27) (Figures S2C and S2D). This result shows that ?-actin mRNA tends to exhibit polarized localization whereas GAPDH mRNA does not. The correlation coefficient between PI and DI describes the distinct localization property of each mRNA species in a population of cells.
% % % 
% % % Full-size image
% % % 
% % %     Figure 2. Comparison of GAPDH and ?-Actin mRNA Distributions in Chicken Embryonic Fibroblasts(A) FISH image of GAPDH mRNA (red), and ?-actin mRNA (green). Scale bar represents 10 ?m.(B) Scatter plot of polarization index for GAPDH mRNA in x axis and ?-actin mRNA in y axis.(C) Scatter plot of dispersion index for GAPDH mRNA in x axis and ?-actin mRNA in y axis. Red dashed lines indicate cells that have the same index for GAPDH and ?-actin mRNA.See also Figure S2.
% % % 
% % % Full-size image
% % % 
% % %     Figure S2. Specificity of Single-Molecule FISH and Quantification of mRNA Distribution, Related to Figure 2(A and B) FISH with MS2_LK51 probe to wild-type MEF (A) and Actb-MBS MEF (B). The intensity range is the same for the two images. Because each ?-actin mRNA in the Actb-MBS MEF contains 12 repeats of the MS2 linker, the bright spots in (B) indicate the presence of MBS-tagged ?-actin mRNA with high specificity. Low fluorescence from nonspecific binding of single probes as shown in (A) is excluded from the analysis by thresholding the intensity value. The specificity of single-molecule FISH method has been also demonstrated previously ( [Femino et al., 1998] and [Zenklusen et al., 2008]). Scale bars, 10 ?m.(C and D) Quantification of ?-actin mRNA and GAPDH mRNA distributions in chicken embryonic fibroblasts. Scatter plot of polarization index in x axis and dispersion index in y axis for ?-actin mRNA (A) and for GAPDH mRNA (B).
% % % 
% % % Finally, we examined the performance of this method to quantify the dynamics of ?-actin mRNA localization in living mammalian cells. Because the correlation between PI and DI for ?-actin mRNA was high, the polarization of mRNA distribution was used as a single metric in live cell analysis. We visualized the endogenous ?-actin mRNA using the cells from the Actb-MBS mouse that contains 24 repeats of MS2 binding site (MBS) cassette in the 3? untranslated region (UTR) of the ?-actin gene (Lionnet et al., 2011). Primary mouse embryonic fibroblasts (MEF) were isolated from the mouse, infected with lentivirus that expresses MS2 capsid protein fused with GFP (MCP-GFP) to allow fluorescent labeling of the ?-actin transcript, and stained with membrane-permeable cytoplasmic dye. A nuclear localization sequence (NLS) was added to the N-terminus of MCP-GFP so that free NLS-MCP-GFPs preferentially localize in the nucleus lowering the background in the cytoplasm (Bertrand et al., 1998). As a negative control, we infected wild-type MEF with NLS-MCP-GFP expressing lentivirus and observed little fluorescence in the cytoplasm (Figure S3A). In the Actb-MBS MEF cells, NLS-MCP-GFP binds to MBS-tagged ?-actin mRNA and the complex is exported out to the cytoplasm (Figure S3B). Using time-lapse imaging, we monitored the localization pattern of ?-actin mRNA as the cell migrated on fibronectin-coated glass surface (Figure 3A and Movie S1). The velocity of the cell centroid was measured at 1-min interval, and defined as the protrusion vector. We found that the mRNA polarization vector (yellow arrows, Figure 3A) and the protrusion vector (red arrows, Figure 3A) were highly correlated in space and in time. In all of the migrating cells that we monitored (n = 11), the autocorrelation curve of protrusion vector decays faster than the one for mRNA polarization vector, indicating that mRNA localization is a slower process than random protrusions. From the mean auto-correlation curves of polarization vector (black curve, Figure 3B) and protrusion vector (red curve, Figure 3B), the half-lives are estimated to be ?16 min for mRNA localization and ?4 min for random protrusions. To quantify the relationship between RNA localization and cell protrusion, we calculated the correlation coefficients of two vectors with varying time lags. Within our time resolution, the mean correlation was the highest at zero lag time (Figure 3C). However, the skewed correlation in the negative time lag suggests that cell protrusion precedes mRNA localization in overall cell movement. It has been shown that mRNA localization is not necessary for random protrusions, but required for directed migration (Shestakova et al., 2001). We found that the net migration distance in 2 hr is highly correlated with the mean polarization index of the cell (r = 0.75) (Figure 3D). This result supports that localization of ?-actin mRNA has a physiological role in directed cell migration.
% % % 
% % % Full-size image
% % % 
% % %     Figure S3. Comparison of NLS-MCP-GFP Fluorescence in Wild-Type MEF and Actb-MBS MEF, Related to Figure 3(A and B) Fluorescence images of NLS-MCP-GFP expressed in wild-type MEF (A) and Actb-MBS MEF (B). The intensity range is the same for the two images. In the absence of MBS, NLS-MCP-GFP localizes to the nucleus as shown in (A). In Actb-MBS MEF, NLS-MCP-GFP bound to MBS-tagged ?-actin mRNA appears in the cytoplasm. Scale bars, 10 ?m.
% % % 
% % % Full-size image
% % % 
% % %     Figure 3. Localization of ?-Actin mRNA in a Migrating Cell(A) Time-lapse images of a mouse embryonic fibroblast. The color map shows the fluorescence intensity of MCP-GFP labeling endogenous ?-actin transcripts divided by the intensity of red fluorescent cytoplasmic dye. Yellow arrows show the polarization vector of mRNA distribution, and red arrows show the protrusion vector of the cell.(B) Mean auto-correlation curves of the polarization vector (black) and the protrusion vector (red) (n = 11 cells). The shaded areas indicate 95% confidence intervals.(C) Mean cross-correlation curve of the polarization vector and the protrusion vector (black curve) and 95% confidence interval (gray area) for n = 11 cells.(D) Migration distance as a function of the time-average of the polarization index over 2 hr. The correlation coefficient between the polarization index and the migration distance is 0.75.See also Figure S3 and Movie S1.
% % % 
% % % In summary, the combination of polarization and dispersion of mRNA distribution enables an effective assessment of mRNA localization. This approach will allow us to quantify subtle changes in RNA distribution by gene deletion or mutation, and to compare the localization characteristics of different transcripts. Moreover, a quantitative analysis of a single cell in time-lapse images will facilitate studies on the kinetics and the physiological role of mRNA localization.
% % % 
% % % We have demonstrated the utility of the new method in budding yeast and primary chicken and mouse embryonic fibroblasts. There is a requirement for this analysis that the entire cell area is imaged while maintaining sufficient sensitivity for detecting RNA signals. In order to apply this method for larger cells such as neurons( [Bassell et al., 1998] and [Lyles et al., 2006]), Drosophila embryos (Lécuyer et al., 2007), and syncytial muscle cells ( [Kislauskis et al., 1993] and [Sigrist et al., 2000]), it may be necessary to image multiple fields in three dimensions and stitch them together to reconstruct the whole specimen (Preibisch et al., 2009). In case of further application to the cells in vivo in the tissue environment, it will be useful to outline the cells by cell-surface markers or plasma membrane stains because the surrounding tissues may obscure the regions of the cell for analysis. The simple unbiasedanalysis of intracellular localization presented in this work may considerably enrich studies on the local regulation and function of many mRNA species.
% % % Experimental Procedures
% % % FISH for Yeast
% % % 
% % % Yeast cells were grown in rich media until they reach early log phase with OD600 ?0.5 (Table S1). Cells were fixed by adding 8 ml of 32% paraformaldehyde to 42 ml of culture for 45 min at room temperature. The rest of the spheroplasting and in situ hybridization followed the procedure described previously ( [Trcek et al., 2012] and [Zenklusen et al., 2008]). Briefly, we synthesized four different single stranded DNA probes each of which was 50-nucleotide-long and labeled with four Cy3 dyes (Table S2). The labeling efficiency of each probe was above 90%, indicating efficient coupling of fluorescent dye with each probe. Cells in G2 and mitosis were selected using morphological markers and analyzed for ASH1mRNA localization (Trcek et al., 2011).
% % % FISH for Chicken Embryonic Fibroblasts
% % % 
% % % Primary chicken embryonic fibroblast (CEF) cells were obtained from Charles River Laboratories (Wilmington, MA). Cells were plated on 10 cm culture dishes and grown at 37°C in MEM supplemented with 10% FCS in an atmosphere of 95% air/5% CO2 for 1–2 days. Cells were trypsinized and seeded on coverslips at a density of 1 × 105 cells/mL. After overnight incubation, cells were washed with DPBS and fixed in 4% paraformaldehyde in PBS for 20 min. The coverslips were stored in 70% ethanol at 4°C for a few days and processed for fluorescence in situ hybridization as described in Singer lab protocols (http://singerlab.org/protocols/). Slides were imaged using an Olympus BX-61 microscope equipped with an X-cite 120 PC lamp (EXFO), a UPlanApo 100× 1.35 NA oil immersion objective (Olympus) and a CoolSNAP HQ CCD camera (Photometrics). We used Chroma filter set 31000 (DAPI), 41007a (Cy3), and 41008 (Cy5). Cells were imaged with 0.2 ?m Z steps in each channel using MetaMorph software (Molecular Devices).
% % % Live Cell Imaging of Mouse Embryonic Fibroblasts
% % % 
% % % Primary MEFs were cultured from 14-day-old embryos isolated from a pregnant female of Actb-MBS mouse (Lionnet et al., 2011). The head and dark cardiac tissue was removed from the embryo and the rest of the body was digested with Trypsin EDTA for 20 min. After adding media (DMEM, 10% FBS, 1% pen/strep), we selected fibroblasts by plating the cells on 10 cm culture dishes for 1 hr and washing off the unattached cells with fresh media. The next day, cells were plated on fibronectin-coated MatTek dishes (MatTek), infected with lentivirus, and incubated for 48 hr to express NLS-MCP-GFP. We stained the cells with 5 ?M CellTracker Orange CMRA (Invitrogen) and replaced the culture media with L-15 media containing 10% FBS, 1% pen/strep, and 1% oxyrase (Oxyrase) prior to the experiment. Time-lapse images were taken on an Olympus IX-71 inverted microscope equipped with a UApo/340 40× 1.35 NA oil immersion objective (Olympus), an MS-2000 XYZ automated stage (ASI) and an iXon electron-multiplying charge-coupled device (EMCCD) camera (Andor). The temperature was kept at 37°C with 60% humidity in an environmental chamber (Precision Plastics). The excitation sources were 488 nm line from an argon ion laser (Melles Griot) and a 561 nm diode-pumped solid state laser (Cobolt). Multiple fields of images were acquired every 1 min using MetaMorph (Molecular Devices).
% % % Image Analysis
% % % 
% % % Image segmentation and quantification were performed using custom designed MATLAB programs. If the cell was larger than the field of view, the 2D/3D Stitching Plugin available through Fiji was used to create a tiled image (Preibisch et al., 2009). To analyze the fixed-cell data, cells were segmented using the DIC or auto-fluorescence image. A binary mask was generated and the centroid of the cell was calculated by the mean value of x- and y-coordinate of the pixels within the mask region. For single-molecule FISH images, a maximum intensity projection of Z stacks was used to quantifymRNA distribution. In budding yeast images, each fluorescent spot was fit with a two-dimensional Gaussian function, yielding the amplitude and the position of the fluorescent particle. From the histogram of fluorescence intensity, the majority of the detected spots were identified as single probes bound nonspecifically. The average intensity of single probes was used to calculate the number of probes bound in each fluorescent spot inside a cell. Only the bright spots binding more than four probes were selected to analyze the number and the position of the target mRNA. The mean and the variance of the x- and y-coordinates of the RNA molecules except for the ones at the transcription sites were used to calculate the polarization and dispersion indices. For FISH images of CEFs, the background was determined from the median intensity value within each cell boundary. After background subtraction, the intensity-weighted centroid and the second moment were calculated.
% % % 
% % % To analyze the live-cell data, masks were generated that define the cell and the nucleus from the CellTracker dye image and the NLS-MCP-GFP image, respectively. The center of the cell was identified by the intensity-weighted centroid of the cytoplasmic dye image. After background subtraction, the NLS-MCP-GFP image was divided by the CellTracker image to obtain RNA distribution image normalized by the cytoplasmic volume. The center of RNA was determined by the weighted centroid of the resulting image excluding the nuclear region. For each time point, the RNA polarization vector was calculated as the vector pointing from the center of the cell to the center of RNA.
% % % 
% % % The protrusion vector was defined as the instantaneous velocity of the cell calculated every 1 min. The net migration distance was determined to be the distance between the beginning and ending points of the cell trajectory. Auto-correlation and cross-correlation of any two sets of vectors were computed by the dot product of the vectors with varied time lags using a similar method demonstrated by Weiger et al. (2010). The half-lives of mRNA polarization and cell protrusion were determined at the correlation coefficient of 0.5 from the autocorrelation curves.
% % % Acknowledgments
% % % 
% % % Microscopy equipment for the live cell imaging experiments was provided by the Gruss Lipper Biophotonics Center. This work was supported by NIH GM57071, GM84364 and GM86217 to R.H.S. H.Y.P. was supported by National Research Service Awards F32-GM087122.
% % % 
% % % Supplemental Information
% % % 
% % %     Table S1. Yeast Strains, Related to Figure 1.  
% % %     Download File (53K)
% % %     Help with PDF files
% % % 
% % %     Table S2. FISH Probes Used in This Study, Related to Figures 1 and 2.  
% % %     Download File (66K)
% % %     Help with PDF files
% % % 
% % %     Movie S1. Localization of ?-actin mRNA during Cell Migration, Related to Figure 3.  Distribution of ?-actin mRNA in a migrating cell in real time. The color map shows the fluorescence intensity of NLS-MCP-GFP divided by the intensity of cytoplasmic dye, representing the concentration of ?-actin mRNA. The time-lapse images were recorded for 3 hr. The localization pattern of ?-actin mRNA is highly dynamic and correlated with the cell protrusion and migration.
% % %     Download Video (3244K)
% % %     Help with AVI files
% % %     
% % %     
    
    
    
    
end