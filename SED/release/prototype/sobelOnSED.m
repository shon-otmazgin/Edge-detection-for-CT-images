function [ sobelOnSEDseg ] = createSEDsegmentation( sedSeg ,folder )
%% evaluate Sobel algoritem on the SED slides
% the function getting the SED seg matrix and dest folder.
% for each slide in the matrix it calculate the Sobel segmentation and return the new sobel matrix.

    workingDir = folder;
    workingDir = fullfile(workingDir,'SED');
    imageNames = dir(fullfile(workingDir,'*.png'));
    imageNames = {imageNames.name}';
    a = size(imageNames);
    length = a(1);
    
    workingDir = folder;
    workingDir = fullfile(workingDir,'SobelOnSED');
    mkdir(workingDir);   %create the directory

    for i=1:length
      
        I = sedSeg(:,:,i);
        E =  edge(I, 'Sobel');
        E = im2single(E);
        sedSeg(:,:,i) = E;
        filename = [sprintf('%03d',i) '.png'];
        fullname = fullfile(workingDir, filename);
        imwrite(E ,fullname);  
    end
    sobelOnSEDseg = sedSeg;
    
end