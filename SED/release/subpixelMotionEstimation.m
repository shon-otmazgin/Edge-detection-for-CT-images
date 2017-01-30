vidReader = VideoReader('prototype\slides\results\SED\SED.avi','CurrentTime',11);

opts.BlockSize   = 5;
opts.SearchLimit = 10;

workingDir = 'prototype\slides\results\SED';
resultSobelDir = 'prototype\slides\results\Subpixel Motion Estimation';
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

figure;
for i=1:length - 1
    frameGrayLast = imread(fullfile(workingDir,imageNames{i}));
    frameGray = imread(fullfile(workingDir,imageNames{i + 1}));
    [MVx, MVy] = Bidirectional_ME(im2double(frameGray), im2double(frameGrayLast)...
                , opts);

    imagesc(MVx.^2);
    set(gca,'XTick',[]); % Remove the ticks in the x axis!
    set(gca,'YTick',[]); % Remove the ticks in the y axis
    set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
    filename = [sprintf('%03d',i) '.jpg'];
    fullname = fullfile(resultSobelDir, filename);
    saveas(gcf,fullname); 
end

%% output a video file.

tic;
videoName = 'Subpixel Motion Estimation.avi';
workingDir = 'prototype\slides\results\Subpixel Motion Estimation';
createVideo(workingDir , videoName);
toc

