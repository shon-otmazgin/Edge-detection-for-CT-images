function [ slide , startFrame, endFrame ] = getOneSlideGroundTruth( scan_name, slideNumber , folder )
% the function is geting liver segmentation name and serach it on MATLAB paths.
% once it found it start to search for the begging og the segmantation
% save the wanted ground truth segmentation and return it.

    workingDir = folder;
    workingDir = fullfile(workingDir,'GroundTruth');
    mkdir(workingDir);   %create the directory
    
    volumeName = strsplit(scan_name,'.');
    segmentationName = strcat(volumeName(1),'_58','.nii.gz');
    segmentationName = strjoin(segmentationName);
    
    CTsegmentation = load_untouch_nii_gzip(segmentationName);
    segMat = CTsegmentation.img;
     z = size(segMat,3);
    startFrame = 0;
    endFrame = 0;
    for ii=1:z
        segmentation = segMat(:,:,ii);
        v = any(segmentation);
        if any(v) ~= 0 && startFrame == 0
            startFrame = ii;
        end
        if any(v) == 0 && startFrame ~= 0
            endFrame = ii - 1;
            break;
        end
    end
    
    segmentation = segMat(:,:,slideNumber);
    
    [~, threshold] = edge(segmentation, 'sobel');
    fudgeFactor = .5;
    slide = edge(segmentation,'sobel', threshold * fudgeFactor); 
    
    filename = [sprintf('%03d',i) '.png'];
    fullname = fullfile(workingDir, filename);
    imwrite(slide ,fullname);
    
   
end

