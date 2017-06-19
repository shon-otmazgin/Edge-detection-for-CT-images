function [ ] = createCTslides(volume,scan_name, folder)
%   the function is getting CT scan name and serach it on MATLAB paths.
%   once the scan it found, the function create slides from the scan and
%   save the slides on the curent folder under 'slides' folder.
    
    scanMat = volume;
    z = size(scanMat,3);
    
    workingDir = folder;
    
    itemFolderName = strsplit(scan_name,'.');
    workingDir = fullfile(workingDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
    workingDir = fullfile(workingDir,'Volume');
    mkdir(workingDir);   %create the directory
   
    for i=1:z
        filename = [sprintf('%03d',i) '.png'];
        fullname = fullfile(workingDir, filename);
        imwrite(scanMat(:,:,i) ,fullname);
    end
end

