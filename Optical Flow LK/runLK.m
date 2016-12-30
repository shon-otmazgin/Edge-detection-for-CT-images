% Video creator
tic;

workingDir = 'test\edges';
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

outputVideo = VideoWriter(fullfile(workingDir,'58_edges.avi'));
outputVideo.FrameRate = 5;  % Default 30
open(outputVideo);

for ii = 200:length
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img);
end
close(outputVideo)

toc