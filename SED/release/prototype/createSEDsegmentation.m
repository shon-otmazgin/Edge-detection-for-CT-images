function [ sedSeg ] = createSEDsegmentation( volume ,folder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    %% set opts for training (see edgesTrain.m)
    opts=edgesTrain();                % default options (good settings)
    opts.modelDir='models/';          % model will be in models/forest
    opts.modelFnm='modelBsds';        % model name
    opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
    opts.useParfor=0;                 % parallelize if sufficient memory

    %% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
    model=edgesTrain(opts); % will load model if already trained

    %% set detection parameters (can set after training)
    model.opts.multiscale=0;          % for top accuracy set multiscale=1
    model.opts.sharpen=2;             % for top speed set sharpen=0
    model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
    model.opts.nThreads=4;            % max number threads for evaluation
    model.opts.nms=0;                 % set to true to enable nms

    %% evaluate edge detector on BSDS500 (see edgesEval.m)
    if(0), edgesEval( model, 'show',1, 'name','' ); end

    %% evaluate SED forest algoritem on the slides
 
    workingDir = folder;
    workingDir = fullfile(workingDir,'Volume');
    imageNames = dir(fullfile(workingDir,'*.png'));
    imageNames = {imageNames.name}';
    a = size(imageNames);
    length = a(1);
    
    workingDir = folder;
    workingDir = fullfile(workingDir,'SED');
    mkdir(workingDir);   %create the directory
    
   
    for i=1:length
      
        I = volume(:,:,i);
        
        %SED input is unit8 hxwx3 image
        I = im2uint8(I);  
        I = cat(3, I, I, I);
        
        E=edgesDetect(I,model);

        volume(:,:,i) = E;
        filename = [sprintf('%03d',i) '.png'];
        fullname = fullfile(workingDir, filename);
        imwrite(E ,fullname);  
    end
    sedSeg = volume;
 end

