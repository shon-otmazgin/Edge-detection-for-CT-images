function [ ] = createCTslides( scan_name, folder)
%UNTITLED Summary of this function goes here
%   the function is getting CT scan name and serach it on 'Volumes' folder.
%   once the scan it found, the function create slides from the scan and
%   save the slides on the curent folder under 'edges' folder.
    CTscan = load_untouch_nii_gzip(scan_name);

    scanMat = CTscan.img;
    z = size(scanMat,3);
    workingDir = folder;
    for i=1:z
        filename = [sprintf('%03d',i) '.jpg'];
        fullname = fullfile(workingDir, filename);
        imwrite(scanMat(:,:,i) ,fullname);
    end
end

