% Demo for Structured Edge Detector (please see readme.txt first).

%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
tic, model=edgesTrain(opts); toc; % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=0;                 % set to true to enable nms

%% evaluate edge detector on BSDS500 (see edgesEval.m)
if(0), edgesEval( model, 'show',1, 'name','' ); end

%% detect edge and visualize results
I = imread('peppers.png');
tic, E=edgesDetect(I,model); toc
figure(1); im(I); figure(2); im(1-E);

% create the EDGES.
%{
workingDir = 'test';
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);
tic;
for ii = 1:length
   img = imread(fullfile(workingDir,imageNames{ii}));
   E=edgesDetect(img,model);
   filename1 = [sprintf('%03d',ii) '.png'];
   fullname1 = fullfile(workingDir,'edges\1',filename1);
   filename2 = [sprintf('%03d',ii) '.png'];
   fullname2 = fullfile(workingDir,'edges\2',filename2);
   imwrite(1-E,fullname1);    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   imwrite(E,fullname2);    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
end
toc
%}

%create Video edges from the images
%{
workingDir = 'test';
imageNames1 = dir(fullfile(workingDir,'edges/1','*.png'));
imageNames1 = {imageNames1.name}';

imageNames2 = dir(fullfile(workingDir,'edges/2','*.png'));
imageNames2 = {imageNames2.name}';

a = size(imageNames1);
length = a(1);

outputVideo1 = VideoWriter(fullfile(workingDir,'1.avi'));
outputVideo1.FrameRate = 5;
open(outputVideo1);
outputVideo2 = VideoWriter(fullfile(workingDir,'2.avi'));
outputVideo2.FrameRate = 5;
open(outputVideo2);
tic;
for ii = 1:length
   img = imread(fullfile(workingDir,'edges\1',imageNames1{ii}));
   writeVideo(outputVideo1,img);
   
   img = imread(fullfile(workingDir,'edges\2',imageNames2{ii}));
   writeVideo(outputVideo2,img);
end
toc
close(outputVideo1)
close(outputVideo2)

%}