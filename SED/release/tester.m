%{
CTscan = load_untouch_nii_gzip('10000014_1_CT_wb.nii.gz');
CTsegmentation = load_untouch_nii_gzip('10000014_1_CT_wb_58_4.nii.gz');
bsrImg = load('BSR\BSDS500\data\groundTruth\train\2092.mat');
%}
%%
workingDir = 'CTtrain\data';

scanMat = CTscan.img;
segMat = CTsegmentation.img;

z = size(scanMat,3);
for i=1:z
    segmentation = segMat(:,:,i);
    v = any(segmentation);
    if any(v) ~= 0
        groundTruth = cell(1 , 7); 
        [~, threshold] = edge(segmentation, 'sobel');
        fudgeFactor = .5;
        boundaries = edge(segmentation,'sobel', threshold * fudgeFactor); 
        for j=1:7
            groundTruth{1, j}.Segmentation = segmentation;
            groundTruth{1, j}.Boundaries = boundaries;
        end
        img.groundTruth = groundTruth;
        filename = [sprintf('%03d',i) '.mat'];
        fullname = fullfile(workingDir,'groundTruth\train', filename);
        save(fullname,'img');
        filename = [sprintf('%03d',i) '.jpg'];
        fullname = fullfile(workingDir,'images\train', filename);
        imwrite(scanMat(:,:,i) ,fullname);
    end
end