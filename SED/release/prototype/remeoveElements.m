clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 15;

tic;

workingDir = 'prototype\slides\results\SED + Sobel';
resultDir = 'prototype\slides\results\detecting cells v3';
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

startFrame = 333; endFrame = 409; 

for i=startFrame:endFrame

     % Read in a demo image.
    grayImage = imread(fullfile(workingDir,imageNames{i}));
    % Get the dimensions of the image.
    % numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
    [rows, columns, numberOfColorChannels] = size(grayImage);
    if numberOfColorChannels > 1
        % It's not really gray scale like we expected - it's color.
        % Use weighted sum of ALL channels to create a gray scale image.
        grayImage = rgb2gray(grayImage);
        % ALTERNATE METHOD: Convert it to gray scale by taking only the green channel,
        % which in a typical snapshot will be the least noisy channel.
        % grayImage = grayImage(:, :, 2); % Take green channel.
    end

    % Binarize the image by thresholding.
    binaryImage = grayImage < 128;

    % Do a moprhological closing on it
    se = strel('disk', 3, 0);
    binaryImage = imclose(binaryImage, se);

    % Extract the largest blob only.
    binaryImage = bwareafilt(binaryImage, 1);

    
  
    filename = [sprintf('%03d',i) '.png'];
    fullname = fullfile(resultDir, filename);
    imwrite(binaryImage ,fullname);
end
toc

%%

tic;
videoName = 'deleting objects.avi';
workingDir = 'prototype\slides\results\detecting cells v3';
createVideo(workingDir , videoName);
toc

