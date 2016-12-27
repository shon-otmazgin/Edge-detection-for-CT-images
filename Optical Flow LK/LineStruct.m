classdef LineStruct
    %LINESTRUCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P1,P2;
        colors;
    end
    methods (Static)
        
        function lineStrctCell = getLineStructCell(overlay,optionalColor,prevLineStrctCell)
            %contourIm - either prior or segmentation
            if size(overlay,4)>1
                boundryIm = sum(overlay,4)>0;
            else
                boundryIm = Utils.getBoundries(overlay);
            end
            
            
            
            
            depth = size(overlay,3);
            if exist('prevLineStrctCell','var')
                lineStrctCell = prevLineStrctCell;
            else
                lineStrctCell = cell(depth,1);
            end
            
            for z = 1:depth
                currentBoundry = boundryIm(:,:,z);
                currentOverlay = squeeze(overlay(:,:,z,:));
                if sum(currentBoundry(:))==0
                    continue;
                end
                if exist('prevLineStrctCell','var') && ~isempty(prevLineStrctCell{z})
                    lineStrctCell{z} = prevLineStrctCell{z};
                else
                    lineStrctCell{z} = LineStruct;
                end
                
                
                CC = bwconncomp(currentBoundry);
                
                for rr = 1:length(CC.PixelIdxList)
                    [yInit,xInit] = ind2sub(size(currentBoundry),CC.PixelIdxList{rr}(1));
                    [yx] = bwtraceboundary(currentBoundry,[yInit xInit],'W',8);
                    origBoundry = [yx(:,1),yx(:,2)];
                    
                    for tt=1:size(origBoundry,1)-1
                        p1 = origBoundry(tt,:);
                        p2 = origBoundry(tt+1,:);
                        %if ~exist('closestPoint1','var')
                        %    closestPoint1 = origBoundry(fMinDist(p1,origBoundry),:);
                        %end
                        %closestPoint2 = origBoundry(fMinDist(p2,origBoundry),:);
                        if exist('optionalColor','var') && ~isempty(optionalColor)
                            color = optionalColor;
                        else
                            c1 = currentOverlay(p1(1),p1(2),:);
                            c2 = currentOverlay(p2(1),p2(2),:);
                            color = squeeze((c1+c2)/2)';
                        end
                        lineStrctCell{z}.P1 = [lineStrctCell{z}.P1 ;p1];
                        lineStrctCell{z}.P2 = [lineStrctCell{z}.P2 ;p2];
                        lineStrctCell{z}.colors = [lineStrctCell{z}.colors ;color];
                        p1 = p2;
                        %closestPoint1 = closestPoint2;
                    end
                end
            end
            
        end
        function out = getContourLineOverlay(contourIm,im,optionalNewColor)
            %contourIm - either prior or segmentation
            boundryIm = Utils.getBoundries(contourIm);
            if length(unique(contourIm(:))) == 2 %case 1 - segmentation
                overlay = Utils.displaySegmentation(contourIm,boundryIm,optionalNewColor);
            else %case 2 - priorResult
                [~,overlay] = Utils.displayCertaintyUncertainty2_3D(im,contourIm,boundryIm);
            end
            
            nLines = 50; %TODO shouldnt be constant
            
            fNorm = @(p1,P) sum((repmat(p1,size(P,1),1)-P).^2,2);
            fMinDist = @(p1,P)find(fNorm(p1,P)==min(fNorm(p1,P)),1,'first');
            %first step - generate RGB image from contourIm
            out = zeros(size(im));
            for z = 1:size(im,3)
                currentSeg = boundryIm(:,:,z);
                currentIm = im(:,:,z);
                currentOverlay = squeeze(overlay(:,:,z,:));
                
                [Y,X] = ind2sub(size(currentSeg),find(currentSeg));
                [yx] = bwtraceboundary(currentSeg,[Y(1) X(1)],'W',8);
                simplifiedBoundry = simplifyPoly(yx, nLines);
                [Y,X] = ind2sub(size(currentSeg),find(currentSeg));
                origBoundry = [Y,X];
                
                
                figure,imshow(currentIm,[]); hold on;
                for tt=1:size(simplifiedBoundry,1)-1
                    p1 = simplifiedBoundry(tt,:);
                    p2 = simplifiedBoundry(tt+1,:);
                    if ~exist('closestPoint1','var')
                        closestPoint1 = origBoundry(fMinDist(p1,origBoundry),:);
                    end
                    closestPoint2 = origBoundry(fMinDist(p2,origBoundry),:);
                    c1 = currentOverlay(closestPoint1(1),closestPoint1(2),:);
                    c2 = currentOverlay(closestPoint2(1),closestPoint2(2),:);
                    color = squeeze((c1+c2)/2);
                    
                    plot([p1(2),p2(2)],[p1(1),p2(1)],'Color',color);
                    
                    p1 = p2;
                    closestPoint1 = closestPoint2;
                end
            end
            
        end
    end
    methods
        
        
        function res = isEmpty(lineStrctObj)
            res = isempty(lineStrctObj.P1);
        end
        
        function plotLines(lineStrctObj)
            hold on;
            for t=1:size(lineStrctObj.P1,1)
                x1 = lineStrctObj.P1(t,2);
                x2 = lineStrctObj.P2(t,2);
                y1 = lineStrctObj.P1(t,1);
                y2 = lineStrctObj.P2(t,1);
                plot([x1,x2],[y1,y2],'Color',lineStrctObj.colors(t,:));
            end
            hold off;
        end
    end
    
end

