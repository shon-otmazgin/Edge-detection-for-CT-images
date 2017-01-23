%% evaluate Sobel algoritem on the SED slides

tic
workingDir = 'prototype\slides\results\SED';
resultSobelDir = 'prototype\slides\results\SED + Sobel';
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

for i=1:length
    I = imread(fullfile(workingDir,imageNames{i}));
    E =  edge(I, 'Sobel');
    filename = [sprintf('%03d',i) '.jpg'];
    fullname = fullfile(resultSobelDir, filename);
    imwrite(1-E ,fullname);
end
toc;

%% output a video file.

tic;
videoName = 'SED + sobel.avi';
workingDir = 'prototype\slides\results\SED + Sobel';
createVideo(workingDir , videoName);
toc