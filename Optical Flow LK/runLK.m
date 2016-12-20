tic;
niiStruct = load_untouch_nii_gzip('10000019_1_CT_wb.nii.gz');

%{
vidReader = VideoReader('lipVid.avi');

opticFlow = opticalFlowLK('NoiseThreshold',0.009);

while hasFrame(vidReader)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);

    flow = estimateFlow(opticFlow,frameGray);

    imshow(frameRGB)
    hold on
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',10)
    hold off
end
%}
toc