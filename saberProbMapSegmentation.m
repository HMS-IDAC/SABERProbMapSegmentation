function sabreProbMapSegmentation(varargin)
%this function requires input of nuclei stack range. It assumes that every
%stack beyond that to the end is a cytoplasmic stain. Marker controlled
%watershed based on distance transform of nuclei channel is employed to
%separate nuclei clumps.


ip = inputParser;
ip.addParamValue('NucMaskChan',[2 2],@(x)(numel(x) == 2 & all(x > 0 )));  
ip.addParamValue('CytoMaskChan',[3 3],@(x)(numel(x) == 2 & all(x > 0 )));  
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','ring'})));
ip.addParamValue('MedianIntensity',true,@islogical);
ip.addParamValue('percentileIntensity',true,@islogical);
ip.addParamValue('percentile',95,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('saveFig',true,@islogical);
ip.addParamValue('saveMasks',true,@islogical);
ip.addParamValue('applyFFC','none',@(x)(ismember(x,{'both','ffonly','none'})));
ip.addParamValue('useRFNuc',true,@islogical);
ip.addParamValue('nuclearRegion','watershed',@(x)(ismember(x,{'distance','watershed','dilation'})));
ip.addParamValue('segmentCytoplasm','segmentCytoplasm',@(x)(ismember(x,{'segmentCytoplasm','loadMask','ignoreCytoplasm'})));
ip.addParamValue('upsample',true,@islogical); % upsamples by a factor of 2 in X and Y when cytoplasm is very thin
ip.addParamValue('useGPUArray',false,@islogical);
ip.parse(varargin{:});          
p = ip.Results;  

%% initialization
drive='D';
parentPath = [drive ':\LSP\Sinem'];
folders=dir([parentPath filesep '*fromPE*']);

FFCPath = [drive ':\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\Real experiment\All Cyles\FFC'];
switch p.applyFFC
    case 'both'
        dfp = volumeRead([FFCPath filesep 'allDF.tif']);         
        ffp = volumeRead([FFCPath filesep 'allFFP.tif']);
    case 'ffonly'
        ffp = volumeRead([FFCPath filesep 'allFFPOnly.tif']);
    case 'none'
end
                            

%%
    for iFolder = 1:numel(folders)
    testPath = [parentPath filesep folders(iFolder).name];
    if exist([testPath filesep 'analysis'],'dir')~=7
        mkdir([testPath filesep 'analysis'])
        
    end
    analysisPath = [testPath filesep 'analysis'];
    
      
              
                    files = dir([testPath filesep  '*_NucSeg.tif']);
                     for iFile = 1:numel(files)
                        
                        tic
                        fileName = files(iFile).name;
                        [~,name,ext] = fileparts(fileName) ;
                        
                        nucleiPM = volumeRead([testPath filesep fileName]);
                        
                        nucleiStack = [1 size(nucleiPM,3)/4];
                        nucleiMaskChan = p.NucMaskChan;
                        
                        nucleiImage = nucleiPM(:,:,nucleiMaskChan(1):nucleiMaskChan(2));
                        nucleiImage = max(nucleiImage,[],3);
                        nucleiImage = double(nucleiImage);
                        nucleiImage = nucleiImage/65535;%max(nucleiImage(:));
                        
                        prefix = fileName(1:strfind(fileName,'NucSeg')-2);
                        channelList = dir([testPath filesep '*' prefix '*']);
                        
                        % read other channels as individual files
                        finalChannelList = [];
                        for iChan = 1:length(channelList)
                            fName = channelList(iChan).name;
                            if ~isdir(fName) && ~contains(fName,'NucSeg')  && ~contains(fName,'._')
                                finalChannelList{end+1} = fName;
%                             
                            end
                        end
                        cytoChanEnd = numel(finalChannelList);

                        if (p.useRFNuc)
                        %% use class prob maps from UNet
                            nucleiClass = imresize(1-nucleiImage,1);

                            %% markers based on log filter on classProbs 3
                             nuclei = nucleiPM(:,:,1);
                             if p.useGPUArray==1
                              logNuclei=  gather(imgaussfilt3(gpuArray(filterLoG(nucleiClass,0.6)),1.1)); %3.5, 2
                             else
                              logNuclei=  imgaussfilt3(filterLoG(nucleiClass,1),1.2); %3.5, 2
                             end
                             logfgm = imregionalmax(logNuclei);
                             threshold= thresholdOtsu(imresize(nuclei,0.1));
                             logfgm = (imresize(nuclei,1)>threshold).*(logfgm==1);

                          %% contours for watershed
                          
                             contours = imresize(nucleiImage,[size(logfgm,1) size(logfgm,2)]);
                        %% apply watershed transform
                        switch p.nuclearRegion
                            case 'watershed'
                              gradmag2= imimposemin(contours,imresize(logfgm,[size(contours,1) size(contours,2)]));
                              foregroundMask= watershed(gradmag2);
                        
                            case 'dilation'
                              foregroundMask = imdilate(logfgm,strel('square',2));
                              
                            case 'distance'
                                logDist = -bwdist(logfgm);
                                cytograd= imimposemin(logDist,logfgm);
                                
                                foregroundMask=watershed(cytograd);
                                nucleiBlur=imgaussfilt3(imresize(nucleiPM(:,:,1),0.5),2);
                                foreground =nucleiBlur>thresholdOtsu(nucleiBlur);
                                foregroundMask = foregroundMask.*cast(imresize(foreground,[size(foregroundMask,1) size(foregroundMask,2)]),class(foregroundMask));
                                        
                         end

                        %% process mask
                       allNuclei = foregroundMask;
                       stats=regionprops(bwlabel(allNuclei),imresize(nucleiClass,[size(allNuclei,1) size(allNuclei,2)]),'MeanIntensity','Area');
                       idx = find([stats.MeanIntensity] > 0.98 & [stats.Area] < 100 );
                       nucleiMask = ismember(bwlabel(allNuclei),idx);
                       
                                             
                       statsNM=regionprops(nucleiMask,'Area');
                       largestNucleusArea=prctile([statsNM.Area],95);
                       

                        else
                       %%use presaved label masks 
                       nucleiMask = imread([testPath filesep name '_nucleiLM' ext]);
                        end
                        
                        %% cytoplasm segmentation
                        
                        switch p.segmentCytoplasm
                            case 'segmentCytoplasm'
                                cyto =[];
                                for iChan = p.CytoMaskChan(1):p.CytoMaskChan(2)
                                   cyto= cat(3,cyto,imread([testPath filesep finalChannelList{iChan}]));
                                end
                                cyto = max(cyto,[],3);
                                
                                switch p.cytoMethod
                                    case 'RF'
                                  F = pcImageFeatures(imresize(double(cyto)/65335,0.5,'bilinear'),modelCat.sigmas,modelCat.offsets,...
                                      modelCat.osSigma,modelCat.radii,modelCat.cfSigma,modelCat.logSigmas,modelCat.sfSigmas,...
                                      modelCat.ridgeSigmas,modelCat.ridgenangs,modelCat.edgeSigmas,modelCat.edgenangs,modelCat.nhoodEntropy,...
                                      modelCat.nhoodStd);
                                     [imL,catClassProbs] = imClassify(F,modelCat.treeBag,100);
                                      contours = imresize(catClassProbs(:,:,2),2);
        %                             bgm = classProbs(:,:,1)>0.9;
        %                             bgm = imresize(bwmorph(bgm,'shrink',Inf),4);
                                      bgm =imresize(bwmorph( imgaussfilt3(cyto,2)<100,'thin',Inf),2);
                                      cytograd= imimposemin(imresize(contours,[size(nucleiMask,1) size(nucleiMask,2)]),bgm|nucleiMask);
                                      cellMask= watershed(cytograd);

                                    case 'distanceTransform'
%                                          bgLoG=imgaussfilt3(filterLoG(max(cyto(:))-cyto,10),10);
%                                          bgMax = imregionalmax(bgLoG);
%                                          threshold = median(bgLoG(bgMax==1));
%                                          bgm=imresize(uint16(bgMax.*(bgLoG<threshold)),2);
%                                        
                                         if p.upsample
                                           nucleiMask = imerode(imresize(nucleiMask,2,'nearest'),strel('square',2));
                                         end
                                         
                                         nMaskDist =-bwdist(~nucleiMask);
                                         cytograd= imimposemin(imresize(nMaskDist,[size(nucleiMask,1) size(nucleiMask,2)]),nucleiMask);
                                         cellMask=watershed(cytograd);
                                         
                                        
                                         thresholdMask=imgaussfilt3(imtophat(imresize(cyto,0.1),strel('disk',100)),3);
                                         threshold = (mode(thresholdMask(:)) + thresholdOtsu(thresholdMask))/2;
                                         preCellMask = bwareaopen(thresholdMask>threshold,400);

                                         cellMask = cellMask.*cast(imresize(preCellMask,[size(cellMask,1) size(cellMask,2)]),class(cellMask));
                                         cellMask = bwlabel(bwareaopen(cellMask,round(largestNucleusArea)));

                                    case 'contours'

                                        contours = normalize(steerableDetector(im2double(cyto),2,1.5));

                                    case 'ring'

                                end

                                nucleiMask = cast(nucleiMask,class(cellMask)).*cellMask;   

                                %clean this up! Eliminates cytoplasms with
                                %no nuclei
                                test=cast(~ismember(unique(cellMask),unique(nucleiMask)),class(unique(cellMask))).*unique(cellMask);
                                bgCells=find(test>0);
                                for iTest=bgCells'
                                    cellMask(cellMask ==test(iTest))=0;
                                end


                                %% eliminate border cells
                                inCells = imclearborder(cellMask>0);
                                borderCells = cellMask.*cast((cellMask>0)-inCells,class(cellMask));
                                borderIdx = unique(borderCells);
                                nucleiMask_border = ~ismember(nucleiMask,borderIdx);

                                cellMask = bwlabel(inCells);
                                nucleiMask = nucleiMask_border.*cellMask;
                                cytoplasmMask = cellMask - nucleiMask;
                
                            case 'loadMask'
                                nucleiMask = bwlabel(nucleiMask);
                                cytoplasmMask = imread([testPath filesep '_' name '_cytoLM' ext]);
                            case 'ignoreCytoplasm'
                                nucleiMask = bwlabel(nucleiMask);
                                cytoplasmMask = nucleiMask;
                                
                        end
                        %% apply flatfield correction using flatfield & darkfield, just flatfield, or no correction 
%                         switch p.applyFFC
%                         
%                             case 'both'
%                                 for iffp = 1:4
%                                     for i =1:(size(nucleiPM,3)/size(ffp,3))
%                                         FFCI(:,:,(iffp-1)*size(nucleiPM,3)/size(ffp,3)+i) = (nucleiPM(:,:,(iffp-1)*size(nucleiPM,3)/size(ffp,3)+i)-dfp(:,:,iffp))./ffp(:,:,iffp) ...
%                                             *mean(reshape(ffp(:,:,iffp),[1 size(ffp,1)*size(ffp,2)]));
%                                     end
%                                 end
%                                 cytoResized=imresize(FFCI,[size(nucleiMask,1) size(nucleiMask,2)]);  
%                                 
%                             case 'ffonly'
%                                 for iffp = 1:4
%                                     for i =1:(size(nucleiPM,3)/size(ffp,3))
%                                         FFCI(:,:,(iffp-1)*size(nucleiPM,3)/size(ffp,3)+i) = (nucleiPM(:,:,(iffp-1)*size(nucleiPM,3)/size(ffp,3)+i))./ffp(:,:,iffp) ...
%                                             *mean(reshape(ffp(:,:,iffp),[1 size(ffp,1)*size(ffp,2)]));
%                                     end
%                                 end
%                                 cytoResized=imresize(FFCI,[size(nucleiMask,1) size(nucleiMask,2)]);  
%                             
%                             case 'none'
%                                 cytoResized= imresize(I,[size(nucleiMask,1) size(nucleiMask,2)]);                        
%                         end

                        %% measure intensities from regions
                        meanIntNucTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        meanIntCytoTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        medianIntNucTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        medianIntCytoTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        prctileIntNucTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        prctileIntCytoTable = zeros(max(nucleiMask(:)),cytoChanEnd);
                        centroidCellTable = zeros(max(nucleiMask(:)),2);

                        for iChan = 1: cytoChanEnd
                            I= imresize(imread([testPath filesep finalChannelList{iChan}]),[size(nucleiMask,1) size(nucleiMask,2)]); 
                            nucleiStats=regionprops(nucleiMask,I,'MeanIntensity','Centroid','Area','PixelIdxList');
                            cytoStats=regionprops(cytoplasmMask,I,'MeanIntensity','Centroid','Area','PixelIdxList');

                            if p.MedianIntensity ==1
                                for iCell = 1: numel(nucleiStats)
                                    medianIntNucTable(iCell,iChan) = median(I(nucleiStats(iCell).PixelIdxList));
                                    medianIntCytoTable(iCell,iChan) = median(I(cytoStats(iCell).PixelIdxList));
                                end
                            end
                            
                            if p.percentileIntensity ==1
                                for iCell = 1: numel(nucleiStats)
                                    prctileIntNucTable(iCell,iChan) = prctile(I(nucleiStats(iCell).PixelIdxList),p.percentile);
                                    prctileIntCytoTable(iCell,iChan) = prctile(I(cytoStats(iCell).PixelIdxList),p.percentile);
                                end
                            end
                        
                            meanIntNucTable(:,iChan) = [nucleiStats.MeanIntensity]';
                            meanIntCytoTable(:,iChan) = [cytoStats.MeanIntensity]';

                        end

                        meanIntTable = [meanIntNucTable meanIntCytoTable];
                        medianIntTable = [medianIntNucTable medianIntCytoTable];
                        prctileIntTable = [prctileIntNucTable prctileIntCytoTable];
                        areaTable = [cat(1,nucleiStats.Area) cat(1,cytoStats.Area)  ];
                        centroidCellTable = cat(1,nucleiStats.Centroid);       

                        %% write results to txt file
                        if ~isempty(areaTable)
                            variableNucNamesMeanIntensity = {};
                            variableCytoNamesMeanIntensity = {};
                            variableNucNamesMedianIntensity = {};
                            variableCytoNamesMedianIntensity = {};
                            variableNucNamesPrctileIntensity = {};
                            variableCytoNamesPrctileIntensity = {};
                            
                            for ivarName = 1:size(meanIntTable,2)/2
                                variableNucNamesMeanIntensity = cat(2,variableNucNamesMeanIntensity,{['NucleiChannelMeanIntensity' int2str(ivarName)]});
                                variableCytoNamesMeanIntensity = cat(2,variableCytoNamesMeanIntensity,{['CytoplasmChannelMeanIntensity' int2str(ivarName)]});
                                variableNucNamesMedianIntensity = cat(2,variableNucNamesMedianIntensity,{['NucleiChannelMedianIntensity' int2str(ivarName)]});
                                variableCytoNamesMedianIntensity = cat(2,variableCytoNamesMedianIntensity,{['CytoplasmChannelMedianIntensity' int2str(ivarName)]});
                                variableNucNamesPrctileIntensity = cat(2,variableNucNamesPrctileIntensity,{['NucleiChannelPrctileIntensity' int2str(ivarName)]});
                                variableCytoNamesPrctileIntensity = cat(2,variableCytoNamesPrctileIntensity,{['CytoplasmChannelPrctileIntensity' int2str(ivarName)]});
                            end
                            writetable(array2table([meanIntTable medianIntTable prctileIntTable areaTable centroidCellTable],'VariableNames',[variableNucNamesMeanIntensity variableCytoNamesMeanIntensity...
                                variableNucNamesMedianIntensity variableCytoNamesMedianIntensity variableNucNamesPrctileIntensity variableCytoNamesPrctileIntensity...
                                'NucleusArea' 'CytoplasmArea' 'CellPosition_X' 'CellPosition_Y']),...
                                [analysisPath filesep name '_cytoMasked.txt'],'Delimiter','\t')
                        end

                        %% display
                         if p.saveMasks==true
                            tiffwriteimj(nucleiMask,[analysisPath filesep name '_nucleiMask.tif'])
                            if p.segmentCytoplasm
                                tiffwriteimj(cytoplasmMask,[analysisPath filesep name '_cytoplasmMask.tif'])
                            end
                         end
                         
                         if p.saveFig==true
                            imshowpair(nucleiMask>0,imresize(sqrt(double(nucleiPM(:,:,1))),[size(nucleiMask,1) size(nucleiMask,2)],'nearest'))
                            savefig ([analysisPath filesep name '_nucleiMasked.fig' ])
                         end
                         close all
                         disp(['Completed ' fileName])
                         toc 
                    end
           
            disp(['Completed all images'])

     end


end