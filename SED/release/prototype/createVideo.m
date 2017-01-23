function [ ] = Untitled2( folder, videoName )
%UNTITLED2 Summary of this function goes here
%   Video creator
    workingDir = folder;
    imageNames = dir(fullfile(workingDir,'*.jpg'));
    imageNames = {imageNames.name}';
    a = size(imageNames);
    length = a(1);

    outputVideo = VideoWriter(fullfile(workingDir, videoName));
    outputVideo.FrameRate = 5;  % Default 30
    open(outputVideo);

    for ii = 200:length
       img = imread(fullfile(workingDir,imageNames{ii}));
       writeVideo(outputVideo,img);
    end
    close(outputVideo)
end

