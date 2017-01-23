classdef Utils
    %UTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
        %input Files
        INPUT_NII_FILES = {'data/volumes/10000100_1_CTce_ThAb.nii.gz'
            'data/volumes/10000104_1_CTce_ThAb.nii.gz'
            'data/volumes/10000105_1_CTce_ThAb.nii.gz'
            'data/volumes/10000106_1_CTce_ThAb.nii.gz'
            'data/volumes/10000108_1_CTce_ThAb.nii.gz'}
        
        INPUT_MAT_FILES = {'data/10000100_1_CTce_ThAb.mat'
            'data/10000104_1_CTce_ThAb.mat'
            'data/10000105_1_CTce_ThAb.mat'
            'data/10000106_1_CTce_ThAb.mat'
            'data/10000108_1_CTce_ThAb.mat'}
        
        GROUND_TRUTHS =   {'data/10000100_1_CTce_ThAb_480_.nii.gz'
            'data/10000104_1_CTce_ThAb_480_.nii.gz'
            'data/10000105_1_CTce_ThAb_480_4.nii.gz'
            'data/10000106_1_CTce_ThAb_480_.nii.gz'
            'data/10000108_1_CTce_ThAb_480_.nii.gz'}
        
        GROUND_TRUTHS_MAT_FILES =   {'data/10000100_1_CTce_ThAb_480_.mat'
            'data/10000104_1_CTce_ThAb_480_.mat'
            'data/10000105_1_CTce_ThAb_480_4.mat'
            'data/10000106_1_CTce_ThAb_480_.mat'
            'data/10000108_1_CTce_ThAb_480_.mat'}
        
        SEGMENTATION_OUT = {'outputDir/10000100_out.nii.gz'
            'outputDir/10000104_out.nii.gz'
            'outputDir/10000105_1_out.nii.gz'
            'outputDir/10000106_1_out.nii.gz'
            'outputDir/10000108_1_out.nii.gz'}
        
        SEGMENTATION_OUT_BONES = {'outputDir/10000100_outb.nii.gz'
            'outputDir/10000104_outb.nii.gz'
            'outputDir/10000105_1_outb.nii.gz'
            'outputDir/10000106_1_outb.nii.gz'
            'outputDir/10000108_1_outb.nii.gz'}
        
        SEGMENTATION_OUT_MAT_FILES = {'outputDir/10000100_out.mat'
            'outputDir/10000104_out.mat'
            'outputDir/10000105_1_out.mat'
            'outputDir/10000106_1_out.mat'
            'outputDir/10000108_1_out.mat'}
        
        BONES_CC_VECT_PREPROCESSING_FILES = {'preprocessing/bonesCCVect10000100_26c.txt'
            'preprocessing/bonesCCVect10000104_26c.txt'
            'preprocessing/bonesCCVect10000105_26c.txt'
            'preprocessing/bonesCCVect10000106_26c.txt'
            'preprocessing/bonesCCVect10000108_26c.txt'}
        
        %first stage - bones segmentation thresholds
        BONES_DEFAULT_MAX_Z = 1300;
        BONES_DEFAULT_MIN_Z_VALS = 150:5:500;
        RADIUS_DIV_FACTOR = 0.55;
        PATCH_SIZE = 30;
        
        ORGAN_TYPE_BONES = 1;
        ORGAN_TYPE_AORTA = 2;
        
    end
    
    methods(Static)
        
        function imOut = imresize3D(im, newSize)
            [y, x, z]= ndgrid(linspace(1,size(im,1),newSize(1)),...
                linspace(1,size(im,2),newSize(2)),...
                linspace(1,size(im,3),newSize(3)));
            imOut=interp3(im,x,y,z);
            
        end
        function [minY,maxY,minX,maxX,minZ,maxZ] = extractMinVals(im,gap)
            if ~exist('gap','var')
                gap = 0;
            end
            
            [Y,X,Z] = ind2sub(size(im),find(im>0));
            [minY,maxY,minX,maxX,minZ,maxZ] = Utils.extractMinValsXYZ(X,Y,Z,size(im),gap);
        end
        
        function [minY,maxY,minX,maxX,minZ,maxZ] = extractMinValsXYZ(X,Y,Z,imSize,gap)
            if ~exist('gap','var')
                gap = 0;
            end
            minX = max(1,min(X)-gap);
            minY = max(1,min(Y)-gap);
            minZ = max(1,min(Z)-gap);
            maxX = min(imSize(2),max(X)+gap);
            maxY = min(imSize(1),max(Y)+gap);
            maxZ = min(imSize(3),max(Z)+gap);
        end
        function savePreProcessingFilesForBonesStage()
            %this function is for debug mode
            for i=1:length(Utils.INPUT_NII_FILES)
                volumeNii = load_untouch_nii_gzip( Utils.INPUT_NII_FILES{i});
                numOfCCVect = Utils.countNumOfCCForThresholds(volumeNii.img, Utils.BONES_DEFAULT_MAX_Z, Utils.BONES_DEFAULT_MIN_Z_VALS);
                save(Utils.BONES_CC_VECT_PREPROCESSING_FILES{i},'numOfCCVect','-ascii');
            end
            
            
        end
        
        function FV = removeIsoSurfaceVertices(FV, wantedNVertices)
            
            if(wantedNVertices > length(FV.vertices))
                error('illegal input');
            elseif wantedNVertices == length(FV.vertices)
                return;
            end
            while (length(FV.vertices) > wantedNVertices)
                FV =reducepatch(FV,length(FV.faces) - 1);
            end
            if  length(FV.vertices) < wantedNVertices
                error('cunt')
            else
                return;
            end
            
            
            %verticesToRemove = [];
            % if(wantedNVertices > length(FV.Vertices))
            %     error('illegal input');
            % elseif wantedNVertices == length(FV.Vertices)
            %     return;
            % end
            
            
            % while isempty(verticesToRemove) || length(verticesToRemove)~=length(unique(verticesToRemove))
            %     nVertices = length(FV.Vertices);
            %     verticesToRemove = randi(nVertices,nVertices-wantedNVertices,1);
            % end
            
            % Remove the vertex values at the specified index values
            % newVertices = FV.Vertices;
            % removedVertices = newVertices(verticesToRemove,:);
            % newVertices(verticesToRemove,:) = [];
            
            % Find the new index for each of the new vertices
            % [~, newVertexIndex] = ismember(FV.vertices, newVertices, 'rows');
            
            % Find any faces that used the vertices that we removed and remove them
            % newFaces = FV.faces(all(FV.faces ~= removedVertices, 2),:);
            
            % Now update the vertex indices to the new ones
            % newFaces = newVertexIndex(newFaces);
            
            % FV.Vertices = newVertices;
            %FV.faces = newFaces;
        end
        
        function outIm = getCurrentPixelCC(im,pixel)
            %receives an image and a pixel. returns an image with the
            %connected components which contains the pixel.
            outIm = zeros(size(im));
            CC = bwconncomp(im);
            for i=1:length(CC.PixelIdxList)
                if sum(CC.PixelIdxList{i}==sub2ind(size(im),pixel(1),pixel(2))) > 0
                    outIm(CC.PixelIdxList{i}) = 1;
                    break
                end
            end
            
        end
        
        function imOut = changeDim(im,dim)
            
            if dim==1
                imOut = zeros(size(im,3),size(im,2),size(im,1));
                for z=1:size(im,1)
                    imOut(:,:,z) = squeeze(im(z,:,:))';
                end
            elseif dim==2
                imOut = zeros(size(im,1),size(im,3),size(im,2));
                for z=1:size(im,2)
                    imOut(:,:,z) = squeeze(im(:,z,:));
                end
            else
                imOut = im;
            end
        end
        
        function mask = linesMask(imgSize, P0, P1, maskIn)
            if exist('maskIn','var');
                mask = maskIn;
            else
                mask = zeros(imgSize);
            end
            
            for ii=1:size(P0,1)
                mask = Utils.lineMask(size(mask),P0(:,ii),P1(:,ii),mask);
            end
        end
        
        function mag = calcGradientMag(im,seg)

            if ~exist('imgradient')
                persistent warnFlag;
                if isempty(warnFlag)
                    warnFlag = true;
                    warning('imgradient not define, using own implementation');
                end
                mag = Prior.imgradient(im);
            else
                mag = zeros(size(im));
                relevantLayers = find(sum(sum(seg,1),2) > 0 & sum(sum(seg,1),2) > 0)';
                for z=relevantLayers
                    mag(:,:,z) = imgradient(im(:,:,z));
                end
            end
            
        end
        function [curvatureArr, xy] = calcCurvatore(seg, halfWinSize,THRESHOLD_MIN,THRESHOLD_MAX)
            if ~exist('THRESHOLD_MIN','var')
                THRESHOLD_MIN = 0;
            end
            if ~exist('THRESHOLD_MAX','var')
                THRESHOLD_MAX = inf;
            end
            
            if ~exist('halfWinSize','var')
                halfWinSize = 5;
            end
            
            seg = Utils.keepBiggestCC(seg);
            [Y,X] = ind2sub(size(seg),find(seg));
            [xy] = bwtraceboundary(seg,[Y(1) X(1)],'W',8);
            
            %curvature calculation
            curvatureArr = zeros(size(xy,1),1);
            for t=1:halfWinSize
                y = xy(t:halfWinSize:end,1);
                x = xy(t:halfWinSize:end,2);
                dx  = gradient(x);
                ddx = gradient(dx);
                dy  = gradient(y);
                ddy = gradient(dy);
                num   = abs(dx .* ddy - ddx .* dy)  + 0.000001;
                denom = dx .* dx + dy .* dy + 0.000001;
                denom = sqrt(denom);
                denom = denom .* denom .* denom;
                curvature = num ./ denom;
                
                %thresholding & normalizing
                curvature(curvature<THRESHOLD_MIN) = 0;
                curvature(curvature>THRESHOLD_MAX) = THRESHOLD_MAX;
                if(max(curvature) > 0)
                    curvature = curvature / max(curvature);
                end
                
                %denoising and max filtering
                %curvature = medfilt1(curvature,3);
                if sum(isnan(curvature)) > 0
                    error('curvature contains nan');
                end
                curvatureArr(t:halfWinSize:end) = curvature;
                
            end
            
        end
        
        function mask = lineMask(imgSize, p0, p1, maskIn)
            
            if exist('maskIn','var');
                mask = maskIn;
            else
                mask = zeros(imgSize);
            end
            p0 = double(p0);
            p1 = double(p1);
            x1 = p0(1); y1 = p0(2);
            x2 = p1(1); y2 = p1(2);
            
            % Distance (in pixels) between the two endpoints
            nPoints = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2));
            
            % Determine x and y locations along the line
            xvalues = round(linspace(x1, x2, nPoints));
            yvalues = round(linspace(y1, y2, nPoints));
            
            % Replace the relevant values within the mask
            mask(sub2ind(size(mask), yvalues, xvalues)) = 1;
            
        end
        
        function outIm = keepBiggestCC(im)
            %receives an image, returned an image with the biggest
            %connected component
            if sum(im(:))==0
                outIm = im;
                return;
            end
            
            numOfCC = bwconncomp(im);
            rp = regionprops(numOfCC, 'Area', 'PixelIdxList');
            [~,ind] = max([rp.Area]);
            outIm = zeros(size(im));
            outIm(rp(ind).PixelIdxList) = 1;
        end
        
        function circleMean = getCircleMean(im, circle)
            %recieves an image and a circle. returns the mean value inside the
            %circle.
            circleMask = Utils.drawCircles(im, circle);
            im(circleMask==0) = 0;
            circleMean = sum(im(:)) / sum(circleMask(:));
        end
        
        function [circleMask overlayIm] = drawCircles(im, peaks)
            %receives an image and a list of circle, returns a mask which
            %contains the circles and an overlay image.
            
            overlayIm = double(im);
            overlayIm = 255*(overlayIm/max(overlayIm(:)));
            overlayIm = repmat(uint8(overlayIm),[1 1 3]);
            
            circleMask = zeros(size(im));
            width = size(im,2);
            height = size(im,1);
            for i=1:size(peaks,2)
                radius = peaks(3,i);
                centerW = peaks(1,i);
                centerH = peaks(2,i);
                [W,H] = meshgrid(1:width,1:height);
                mask = (sqrt((W-centerW).^2 + (H-centerH).^2) - radius)<1;
                circleMask(mask) = 1;
                tempRed = overlayIm(:,:,1);
                tempRed(mask) = 255;
                overlayIm(:,:,1) = tempRed;
            end
            
        end
        
        function [patchAroundCenter patchTopLeft] = extractCenterPatch(binaryIm, center)
            %receives a binary image and a center and extract a patch
            %from the center of the image.
            
            y = center(1);
            x = center(2);
            patchAroundCenter = zeros(size(Utils.PATCH_SIZE*2+1,Utils.PATCH_SIZE*2+1));
            patchTopLeft = [0 0];
            if x~=0 && y~=0
                patchTopLeft = [y-Utils.PATCH_SIZE,x-Utils.PATCH_SIZE];
                patchBottomRight = [y+Utils.PATCH_SIZE,x+Utils.PATCH_SIZE];
                patchAroundCenter = binaryIm(patchTopLeft(1):patchBottomRight(1),patchTopLeft(2):patchBottomRight(2));
            end
        end
        
        function [centerPixel, patchAroundCenter patchTopLeft] = extractCenter(binaryIm)
            %receice an image and find the center of the object within the image
            contour = imdilate(binaryIm,ones(3,3)) - binaryIm;
            contourInd = find(contour);
            [Y X] = ind2sub(size(binaryIm),contourInd);
            contour(1:min(Y),:) = 1;
            contour(max(Y):end,:) = 1;
            contour(:,1:min(X)) = 1;
            contour(:,max(X):end) = 1;
            
            distFromContour = bwdist(contour);
            [~,centralIndx] = max(distFromContour(:));
            [Cy Cx] = ind2sub(size(binaryIm),centralIndx);
            centerPixel = [Cy, Cx];
            [patchAroundCenter patchTopLeft] = Utils.extractCenterPatch(binaryIm,centerPixel);
        end
        
        function dist = calcCenterDistUniformity(patch)
            %receive a patch. calculates variance of distances from the center, along it's contour.
            center = int16(ceil(size(patch)/2));
            centerMask = zeros(size(patch));
            centerMask(center(1),center(2)) = 1;
            distFromCenter = bwdist(centerMask);
            
            contour = imdilate(patch,ones(3,3)) - patch;
            distAlongContour = distFromCenter(logical(contour));
            
            dist = std(distAlongContour);
        end
        
        function dist = calcDistBetweenContours(patch1,patch2)
            %calculates a distance function between two patches.
            contour1 = imdilate(patch1,ones(3,3)) - patch1;
            contour2 = imdilate(patch2,ones(3,3)) - patch2;
            distFromC1 = bwdist(contour1);
            distC2FromC1 = distFromC1(logical(contour2));
            dist = max(distC2FromC1);
        end
        
        function [peaks,circleMask,overlayIm] = findCircles(im, nCircles)
            %receives an image and a number nCircles. finds nCircle circles
            %in the image.
            if ~exist('nCircles','var')
                nCircles = 7;
            end
            e = edge(im, 'canny');
            radii = 10:1:22;
            h = circle_hough(e, radii, 'same', 'normalise');
            peaks = circle_houghpeaks(h, radii, 'nhoodxy', 3, 'nhoodr', 21, 'npeaks', nCircles);
            
            [circleMask overlayIm] = Utils.drawCircles(im, peaks);
        end
        
        function imOut = imfillSlices(imIn)
            %receives a 3D binary image. for every 2D slice in the image -
            %fills the holes inside it.
            imOut = imIn;
            for z=1:size(imIn,3)
                imOut(:,:,z) = imfill(imIn(:,:,z));
            end
        end
        
        function out = getBoundries(mask, uesBrokenLines)
            if ~exist('uesBrokenLines','var')
                uesBrokenLines = false;
            end
            
            out = zeros(size(mask));
            for z=1:size(mask,3)
                curMask = mask(:,:,z);
                boundryPoints  = bwboundaries(curMask,8);
                
                
                curOut = zeros(size(curMask));
                
                for i=1:length(boundryPoints)
                    if uesBrokenLines
                        boundryPoints{i} = boundryPoints{i}(1:2:end,:);
                    end
                    inds =  sub2ind(size(curMask),boundryPoints{i}(:,1),boundryPoints{i}(:,2));
                    curOut(inds) = 1;
                end
                curOut = logical(curOut);
                out(:,:,z) = curOut;
            end
            out = logical(out);
        end
        
        
        
        function [outOverlay,outColor] = displayCertaintyUncertainty2_3D(im,map,segMask, gThresh,rThresh)
            
            if ~exist('gThresh','var')
                gThresh = 1;
            end
            if ~exist('rThresh','var')
                rThresh = 0;
            end
            
            if max(im(:))>0
                for z=1:size(im,3)
                    im(:,:,z) = mat2gray(im(:,:,z));
                end
            end
            
            outColor = zeros(size(map,1),size(map,2),size(im,3),3);
            outOverlay = zeros(size(map,1),size(map,2),size(im,3),3);
            
            for z=1:size(im,3)
                [overlay,color] = Utils.displayCertaintyUncertainty2(im(:,:,z),map(:,:,z),segMask(:,:,z),gThresh,rThresh);
                outColor(:,:,z,:) = color;
                outOverlay(:,:,z,:) = overlay;
            end
        end
        
        function displaySegmentation2D(im,seg)
            im = mat2gray(im);
            boundries = Utils.getBoundries(seg);
            figure,imshow(im);hold on;
            [Y,X] = find(boundries);
            plot(X,Y,'.');
            
            
        end
        
        
        function [out,valRange,indRange] = optimizeImContrast(im)

            f = figure
            h = histogram(im);
            LOW_PERCENTAGE_THRESH = 0.001*prod(size(im));
            %finds max peak, searches for upper and lower thresholds.
            [~,maxInd] = max(h.Values);
            
            lowPrecentageMask = h.Values<LOW_PERCENTAGE_THRESH;
            
            localMinimaMask = conv(h.Values,[0,0,0,0,0,0,1,0,0,0,0,0,-1],'same')<0 & conv(h.Values,[-1,0,0,0,0,0,1,0,0,0,0,0,0],'same')<0;
            finalMask = lowPrecentageMask | localMinimaMask;
            leftBound = find(finalMask(1:maxInd),1,'last');
            rightBound = maxInd+find(finalMask(maxInd+1:end),1,'first');
            valRange=  [h.BinEdges(leftBound),h.BinEdges(rightBound)];
            indRange = [leftBound,rightBound];
            out = mat2gray(im,valRange);
            close(f);
            
        end
        
        function out = displaySegmentation(im,seg,rgbVal,useBrokenLines)
            
            if ~exist('rgbVal','var')
                rgbVal = [1 0 0];
            end
            if ~exist('useBrokenLines','var')
                useBrokenLines = false;
            end
            
            boundries = Utils.getBoundries(seg,useBrokenLines);
            out = zeros(size(im,1),size(im,2),size(im,3),3);
            
            for z = 1:size(im,3)
                
                
                if ndims(im)==4
                    currentIm = squeeze(im(:,:,z,:));
                    out(:,:,z,1) = currentIm(:,:,1);
                    out(:,:,z,2) = currentIm(:,:,2);
                    out(:,:,z,3) = currentIm(:,:,3);
                else
                    currentIm = mat2gray(im(:,:,z));
                    out(:,:,z,1) = currentIm;
                    out(:,:,z,2) = currentIm;
                    out(:,:,z,3) = currentIm;
                end
                
                
                if sum(sum(seg(:,:,z)))==0
                    continue;
                end
                
                rChannel =  out(:,:,z,1);
                gChannel = out(:,:,z,2);
                bChannel = out(:,:,z,3);
                rChannel(boundries(:,:,z)) = rgbVal(1);
                gChannel(boundries(:,:,z)) = rgbVal(2);
                bChannel(boundries(:,:,z)) = rgbVal(3);
                
                out(:,:,z,1) = rChannel;
                out(:,:,z,2) = gChannel;
                out(:,:,z,3) = bChannel;
            end
        end
        
        function imOut = interpolateIm(im,roi2Fill)
            roi2Fill(1,:) = 0;
            roi2Fill(end,:) = 0;
            roi2Fill(:,1) = 0;
            roi2Fill(:,end) = 0;
            
            [X, Y] = ind2sub(size(roi2Fill),find(~roi2Fill));
            [Xq, Yq] = ind2sub(size(roi2Fill),find(roi2Fill));
            V = im(~roi2Fill);
            
            Vq = interp2(X,Y,V,Xq,Yq);
            im(sub2ind(size(roi2Fill),X,Y)) = Vq;
            imOut = im;
        end
        
        function [outOverlay,outColor] = displayCertaintyUncertainty2(im,map,segMask, gThresh,rThresh)
            if exist('gThresh','var')
                map(map>=gThresh) = 1;
            end
            if exist('rThresh','var')
                map(map<=rThresh) = 0;
            end
            if ~exist('segMask','var')
                segMask = Utils.getBoundries(map>=0);
            end
            
            
            QUANT = 16;
            R = linspace(1,0,QUANT+1); %+1 is for handling 0
            G = linspace(0,1,QUANT+1);
            B = zeros(size(R));
            map = ceil(map*16);
            
            Rc = zeros(size(map));
            Gc = zeros(size(map));
            Bc = zeros(size(map));
            for i=0:QUANT
                Rc(map == i & segMask) = R(i+1);%+1 is for handling 0
                Gc(map == i & segMask) = G(i+1);
                Bc(map == i & segMask) = B(i+1);
            end
            outColor = zeros(size(map,1),size(map,2),3);
            outColor(:,:,1) = Rc;
            outColor(:,:,2) = Gc;
            outColor(:,:,3) = Bc;
            
            %# create a gaussian mask for transparency
            %[r,c,~] = size(outColor);
            %M = imgaussfilt(double(segMask), 2);
            %M = mat2gray(M);
            
            %# show overlayed images
            %I_transp  = outColor;
            %figure, imshow(im, 'XData',[1 c], 'YData',[1 r]), hold on
            %hImg = imshow(I_transp);
            %set(hImg, 'AlphaData',M);
            
            
            outOverlay = double(repmat(im,[1,1,3]));
            out2R = double(im); Out2B=double(im); Out2G=double(im);
            out2R(segMask) = Rc(segMask);
            Out2B(segMask) = Bc(segMask);
            Out2G(segMask) = Gc(segMask);
            outOverlay(:,:,1) = out2R;
            outOverlay(:,:,2) = Out2G;
            outOverlay(:,:,3) = Out2B;
            
            
        end
        
        function out = displayCertaintyUncertainty(im,map,segMask,T)
            certainty = map > T & segMask;
            uncertainty = map < T & segMask;
            
            certainty3D = repmat(certainty,[1 1 3]);
            
            imCertaintyOverlay = Utils.combine2ImColor(im,certainty,[],'b');
            imUncertaintyOverlay = Utils.combine2ImColor(im,uncertainty,[],'r');
            
            out = imCertaintyOverlay;
            out(~certainty3D) = imUncertaintyOverlay(~certainty3D);
            
            %imshow(out);
        end
        
        function outIm = combine2ImColor(im1,im2,mask,color)
            if ~exist('mask','var') || isempty(mask)
                mask = im2 > 0;
            end
            if ~exist('color','var')
                layer = 1;
            elseif isequal(color,'r')
                layer =1;
            elseif isequal(color,'g')
                layer =2;
            else
                layer = 3;
            end
            %for debug porpuses - combine two 3D combine2ImColormasks together.
            im1 = uint8(mat2gray(im1)*255);
            im2 = uint8(mat2gray(im2)*255);
            
            outIm = repmat(im1,[1,1,3]);
            mask3d = repmat(mask,[1,1,3]);
            outIm(mask3d) = 0;
            
            markLayer = im2;
            markLayer(~mask) = im1(~mask);
            outIm(:,:,layer) = markLayer;
            
        end
        
        function rect = getBoundingBox(seg)
            rect = zeros(size(seg));
            for z=1:size(seg,3)
                [Y,X] = ind2sub(size(seg(:,:,z)),find(seg(:,:,z)));
                rect(min(Y):max(Y),min(X):max(X),z) = 1;
            end
            %rect(min(Y):max(Y),min(X):max(X),min(Z):max(Z)) = 1;
        end
        
        function rect = getConvHullBox(seg)
            rect = zeros(size(seg));
            for z=1:size(seg,3)
                currSeg = seg(:,:,z);
                rect(:,:,z) = bwconvhull(currSeg);
                
            end
            %rect(min(Y):max(Y),min(X):max(X),min(Z):max(Z)) = 1;
        end
        
        function [outImCell] = loadNiiFiles(strCell,maskFlag)
            
            if ~exist('maskFlag','var')
                maskFlag = false;
            end
            for ii=1:length(strCell)
                outImCell{ii} = IO.loadFile(strCell{ii});
                if maskFlag
                    outImCell{ii}.img = logical(outImCell{ii}.img);
                end
            end
        end
        
        function [outImCell] = loadImages(strCell,maskFlag)
            
            if ~exist('maskFlag','var')
                maskFlag = false;
            end
            outImCell = cell(size(strCell));
            for ii=1:length(strCell)
                outImCell{ii} = IO.loadFile(strCell{ii});
                if isstruct(outImCell{ii})
                    outImCell{ii} = outImCell{ii}.img;
                end
                
                if maskFlag
                    outImCell{ii} = logical(outImCell{ii});
                end
            end
            
        end
        
        function outIm = combine2Im(im1,im2)
            %for debug porpuses - combine two 3D masks together.
            outIm = zeros(size(im1));
            im1 = im1 / max(im1(:));
            for z=1:size(im1,3)
                outIm(:,:,z) = im1(:,:,z) + 0.15*im2(:,:,z);
            end
            
        end
        
        function outIm = concatImages(imCell, horizontalCat)
            if ~exist('horizontalCat','var')
                horizontalCat = true;
            end
            outIm =  imCell{1};
            for i=2:length(imCell)
                outIm = Utils.concat2Im(outIm,imCell{i},horizontalCat);
            end
            
        end
        
        function outIm = concat2Im(im1,im2, horizontalCat)
            %for debug porpuses - combine two 3D masks together.
            if ~exist('horizontalCat','var')
                horizontalCat = true;
            end
            im1 = im1 / max(im1(:));
            im2 = im2 / max(im2(:));
            [h1, w1, z, t] = size(im1);
            [h2, w2,~,~] = size(im2);
            
            if horizontalCat
                %outIm = zeros([h1,w1+w2,z,t]);
                %for z=1:size(im1,3)
                %    outIm(:,1:w1,z,:) = im1(:,:,z,:);
                %    outIm(:,w1+1:end,z,:) = im2(:,:,z,:);
                %end
                outIm = cat(2,im1,im2);
            else
                %outIm = zeros([h1+h2,w1,z,t]);
                %for z=1:size(im1,3)
                %    outIm(1:h1,:,z,:) = im1(:,:,z,:);
                %    outIm(h1+1:end,:,z,:) = im2(:,:,z,:);
                %end
                outIm = cat(1,im1,im2);
            end
        end
        
        function boundary = getBoundaryMask(seg)
            error('use getboundry');
            boundary = zeros(size(seg));
            
            for z=1:size(seg,3)
                currentSeg = seg(:,:,z);
                pointsStrct = bwboundaries(currentSeg,8,'holes');
                points = pointsStrct{1};
                currentSeg = seg(:,:,z);
                for t=1:length(pointsStrct)
                    ind = sub2ind(size(currentSeg),points(:,1),points(:,2));
                    mask = zeros(size(currentSeg));
                    mask(ind) = 1;
                    boundary(:,:,z) = boundary(:,:,z) | mask;
                    
                end
            end
        end
        
        function imOut = viewXz(imIn, type)
            %for debug porpuses - viewing a 3D image in different axis.
            if type=='x'
                imOut = zeros([size(imIn,1) size(imIn,3) size(imIn,2)]);
                for i=1:size(imIn,2)
                    imOut(:,:,i) = reshape(imIn(:,i,:),size(imOut(:,:,i)));
                end
            elseif type=='y'
                imOut = zeros([size(imIn,2) size(imIn,3) size(imIn,1)]);
                for i=1:size(imIn,1)
                    imOut(:,:,i) = reshape(imIn(i,:,:),size(imOut(:,:,i)));
                end
                imOutRot = zeros([size(imOut,2) size(imOut,1) size(imOut,3)]);
                for i=1:size(imOut,3)
                    imOutRot(:,:,i) = imrotate(imOut(:,:,i),90);
                end
                imOut = imOutRot;
            end
            figure,imshow3D(imOut);
        end
        
        
        function count = jointHistogram(image1, bins1, image2, bins2)
            %JOINTHISTOGRAM Summary of this function goes here
            %   calculates the joint histogram of two images.
            %   image1 and image2 are 3D matrices.
            %   bins1 and bins2 describe the histogram bins.
            %   bins are monotonic increased.
            %   example: [0 50 100 150] -> [0,50), [50,100), [100,150)
            %   output: count is 2D histogram with size
            %   (length(bin1)-1)x(length(bin2)-1)
            
            %intializes count
            count = zeros(length(bins1)-1, length(bins2)-1);
            EPSILON = 0.00001;
            bins2(end) = bins2(end) - EPSILON;
            %precalculated all the masks for image 2 (if there is eonugh memory).
            try
                image1RangeMasks = cell(length(bins1)-1,1);
                for i=1:length(bins1)-1
                    mask = getMaskBetweenRange(image1,bins1(i),bins1(i+1));
                    image1RangeMasks{i} = mask;
                end
                image1RangeMasksInitialized = true;
            catch
                clear image2RangeMasks;
                image1RangeMasksInitialized = false;
            end
            
            %calculates the joint histogram
            for i=1:(length(bins1)-1)
                if image1RangeMasksInitialized
                    im1Mask = image1RangeMasks{i};
                else
                    im1Mask = Utils.getMaskBetweenRange(image1,bins1(i),bins1(i+1));
                end
                
                im2RelevanceMask = image2;
                im2RelevanceMask(im1Mask==0) = bins2(1)-1;
                count(i,:) = histcounts(im2RelevanceMask(:),bins2);
            end
            
        end
        
        
        function mask = getMaskBetweenRange(im, minVal, maxVal)
            %Auxiliary function. gets an image and a grayscale range,
            %return a mask of all the pixels within this range
            mask = (im >= minVal) & (im < maxVal);
        end
        
        function score = mutualInformation(image1,image2)
            %MUTUALINFORMATION Summary of this function goes here
            %   Detailed explanation goes here
            % calculates the mutual informatoin of two input images.
            % MI(im1,im2) = H(im1) + H(im2) - H(im1,im2)
            % when H(A) = -sum_i(p(ai).*log(p(ai)))
            %      H(A,B) = -sum_ij(p(ai,bj)log(p(ai,bj)))
            
            %calculates histograms
            nbins = 140;
            [h1, range1] = histcounts(image1,nbins);
            [h2, range2] = histcounts(image2,nbins);
            h12 = Utils.jointHistogram(image1,range1,image2,range2);
            h12 = h12(:);
            
            %implies the probabilities
            h1 = h1 / sum(h1(:));
            h2 = h2 / sum(h2(:));
            h12 = h12 / sum(h12(:));
            
            %removes zeros
            h1(h1==0) = [];
            h2(h2==0) = [];
            h12(h12==0) = [];
            
            %calculates entropies
            H1 = -sum(h1.*log2(h1));
            H2 = -sum(h2.*log2(h2));
            H12 = -sum(h12.*log2(h12));
            
            %returns final results
            score = H1 + H2 - H12;
            
        end
        
        
        
        
        function [symmAxe topSymm symmMatDisp] = getSymmetyAxe(im3D, minZ, maxZ)
            %receive a 3D binary image with z layers, returns a list of
            %size z, which contains the symmetryAxe in the image.
            SIMMETRICALITY_WEIGHT = 1;
            BOTTOM_TOP_SUM_DIFF_WEIGHT = 1;
            LENGTH_WEIGHT = 1;
            
            allSymmetricalAxis = zeros(maxZ-minZ+1,1);
            for z=minZ:maxZ
                im = im3D(:,:,z);
                sumOfRows = sum(im,2);
                minY = find(sumOfRows,1,'first');
                maxY = find(sumOfRows,1,'last');
                
                WIN_SIZE = min(minY-1,size(im,1)-maxY);
                
                allLengths = zeros(maxY-minY+1,1);
                allSimmVals = zeros(maxY-minY+1,1);
                allBottomTopsDiffs = zeros(maxY-minY+1,1);
                for y = minY:maxY
                    lineLen = find(im(y,:),1,'last')-find(im(y,:),1,'first') + 1;
                    symetricalityRes = sum(sum(abs(im(y-WIN_SIZE:y,:)-im(y+WIN_SIZE:-1:y,:))));
                    allBottomTopsDiff = abs(sum(sum(im(y-WIN_SIZE:y,:))) - sum(sum(im(y:y+WIN_SIZE,:))));
                    allLengths(y-minY+1) = lineLen;
                    allSimmVals(y-minY+1) = 1/symetricalityRes;
                    allBottomTopsDiffs(y-minY+1) = 1/allBottomTopsDiff;
                end
                
                allLengths = allLengths / max(allLengths);
                allSimmVals = allSimmVals / max(allSimmVals);
                allBottomTopsDiffs = allBottomTopsDiffs / max(allBottomTopsDiffs);
                
                finalScoreVect = SIMMETRICALITY_WEIGHT*allSimmVals + LENGTH_WEIGHT*allLengths;
                finalScoreVect = finalScoreVect + BOTTOM_TOP_SUM_DIFF_WEIGHT*allBottomTopsDiffs;
                [~, yOut] = max(finalScoreVect);
                yOut + minY - 1;
                allSymmetricalAxis(z-minZ+1) = yOut + minY - 1;
                
            end
            allSymmetricalAxis = medfilt1(allSymmetricalAxis,5);
            allSymmetricalAxis = allSymmetricalAxis + 20; %adds 10 pixels just in case
            symmMatDisp = im3D;
            topSymm = im3D;
            symmAxe = zeros(size(topSymm,3),1);
            for z=1:length(allSymmetricalAxis)
                currentZ = z + minZ - 1;
                ySymm = allSymmetricalAxis(z);
                symmAxe(currentZ) = ySymm;
                topSymm(ySymm:end,:,currentZ) = 0;
                topSymm(1:ySymm-1,:,currentZ) = 1;
                
                currentWristCpy = symmMatDisp(ySymm:end,:,currentZ);
                currentWristCpy(currentWristCpy==0) = 2;
                symmMatDisp(ySymm:end,:,currentZ) = currentWristCpy;
                
            end
        end
        
    end
    
end

