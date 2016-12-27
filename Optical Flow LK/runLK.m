tic;
%niiStruct = load_untouch_nii_gzip('10000067_1_CT_wb.nii.gz');
ct_img = niiStruct.img;

vidReader = VideoReader('test.avi');

%opticFlow = opticalFlowLK('NoiseThreshold',0.009);
opticFlow = opticalFlowHS;

while hasFrame(vidReader)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);

    flow = estimateFlow(opticFlow,frameGray);

    imshow(frameRGB)
    hold on
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',25)
    hold off
end


% Video creator
%{
workingDir = 'test';
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
open(outputVideo);

for ii = 1:length
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img);
end
close(outputVideo)
%}

toc