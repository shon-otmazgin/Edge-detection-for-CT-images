function [ sobelOnSEDseg ] = createSEDsegmentation( sedSeg ,folder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% evaluate Sobel algoritem on the SED slides

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