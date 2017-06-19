function [ finalSeg ] = traceLiver( segMat ,slide ,stratFrame, endFrame, GTframe, folder )
%the function getting the sobel seg matrix, start frame number segmentation
% end frame numer segmentation, ground truth silde and floder.
% for each slide in the sobel segmentation it calculate the active contour using the slide as a mask

    workingDir = folder;
    workingDir = fullfile(workingDir,'SobelOnSED');
    imageNames = dir(fullfile(workingDir,'*.png'));
    imageNames = {imageNames.name}';
    a = size(imageNames);
    l = a(1);
    
    originalSlide = slide;
    
    workingDir = folder;
    workingDir = fullfile(workingDir,'FinalSegmentation');
    mkdir(workingDir);   %create the directory
    
    for ii=1:stratFrame-1
        segMat(:,:,ii) = 0;
    end
    
    slide = imfill(originalSlide, 'holes');
    segMat(:,:,GTframe) = slide;
    
    for ii=GTframe-1:-1:stratFrame

        %open the next slide
        nextSlide = segMat(:,:,ii);
        sh = nextSlide;
        % create the mask
        mask = imdilate(slide ,strel('disk' , 15));
        
        %outside the mask countor set 0.
        nextSlide(mask == 0 ) = 0;
        
        % replace the values in order to get active contour
        nextSlide(nextSlide == 0 ) = 255;
        nextSlide(nextSlide == 1 ) = 0;
        
        % Active contur
        bw = activecontour(nextSlide,mask,500,'edge');
        
        % return to original values - not nesseary.
        nextSlide(nextSlide == 0 ) = 1;
        nextSlide(nextSlide == 255 ) = 0;
        
        %save the segmentation
        
        [B,L] = bwboundaries(bw,'noholes');
        if length(B) ~= 0
            boundary = B{1};

            z = length(boundary);
            mat = zeros(size(bw));
            for iii=1:z
                mat(int16(boundary(iii,1)), int16(boundary(iii,2)))= 1; 
            end
            bw = imfill(mat, 'holes');    
        else
            bw = 0;
            bw = im2single(bw);
            segMat(:,:,ii) = bw;
            for iiii=ii-1:-1:stratFrame
                segMat(:,:,iiii) = 0;
            end
            break;
        end
        
        bw = im2single(bw);
        segMat(:,:,ii) = bw;
        
        % go to next slide.
        slide = bw;

        %disp(ii)
        filename = [sprintf('%03d',ii) '.png'];
        fullname = fullfile(workingDir, filename);
        imwrite(bw ,fullname); 
        
    end
    
    slide = imfill(originalSlide, 'holes');
    
    for ii=GTframe+1:endFrame

        %open the next slide
        nextSlide = segMat(:,:,ii);
        sh = nextSlide;
        % create the mask
        mask = imdilate(slide ,strel('disk' , 15));
        
        %outside the mask countor set 0.
        nextSlide(mask == 0 ) = 0;
        
        % replace the values in order to get active contour
        nextSlide(nextSlide == 0 ) = 255;
        nextSlide(nextSlide == 1 ) = 0;
        
        % Active contur
        bw = activecontour(nextSlide,mask,500,'edge');
        
        % return to original values - not nesseary.
        nextSlide(nextSlide == 0 ) = 1;
        nextSlide(nextSlide == 255 ) = 0;
        
        
        %save the segmentation
        [B,L] = bwboundaries(bw,'noholes');
        if length(B) ~= 0
            boundary = B{1};

            z = length(boundary);
            mat = zeros(size(bw));
            for iii=1:z
                mat(int16(boundary(iii,1)), int16(boundary(iii,2)))= 1; 
            end
            bw = imfill(mat, 'holes');    
        else
            bw = 0;
            bw = im2single(bw);
            segMat(:,:,ii) = bw;
            for iiii=ii+1:-1:endFrame
                segMat(:,:,iiii) = 0;
            end
            break;
        end  
        
        bw = im2single(bw);
        segMat(:,:,ii) = bw;
        
        % go to next slide.
        slide = bw;
        
        %disp(ii)
        filename = [sprintf('%03d',ii) '.png'];
        fullname = fullfile(workingDir, filename);
        imwrite(bw ,fullname); 
        
    end
    
    z = size(segMat,3);
    for ii=endFrame+1:z
        segMat(:,:,ii) = 0;
    end
    
    finalSeg = segMat;
end