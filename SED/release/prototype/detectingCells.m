tic
workingDir = 'prototype\slides\results\SED';
resultDir = 'prototype\slides\results\detecting cells';
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

for i=361:361%length
    I = imread(fullfile(workingDir,imageNames{i}));
    [~, threshold] = edge(I, 'sobel');
    fudgeFactor = .5;
    BWs = edge(I,'sobel', threshold * fudgeFactor);
    %figure, imshow(BWs), title('binary gradient mask');
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    BWsdil = imdilate(BWs, [se90 se0]);
    %figure, imshow(BWsdil), title('dilated gradient mask');
    BWdfill = imfill(BWsdil, 'holes');
    %figure, imshow(BWdfill);
    title('binary image with filled holes');
    BWnobord = imclearborder(BWdfill, 4);
    %figure, imshow(BWnobord), title('cleared border image');
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    %figure, imshow(BWfinal), title('segmented image');
    
    imshow(BWfinal)
  
    str = 'Click to select initial contour location. Double-click to confirm and proceed.';
    title(str,'Color','b','FontSize',12);
    disp(sprintf('\nNote: Click close to object boundaries for more accurate result.'))
    
    mask = roipoly;
  
    figure, imshow(mask)
    title('Initial MASK');
    
    maxIterations = 200; 
    bw = activecontour(BWfinal, mask, maxIterations, 'Chan-Vese');

    % Display segmented image
    figure, imshow(bw)
    title('Segmented Image');
    

    
   
    
    %{
    
    %}
    
 
    

   
     
    filename = [sprintf('%03d',i) '.jpg'];
    fullname = fullfile(resultDir, filename);
    imwrite(bw ,fullname);
end
toc;


%%

tic;
videoName = 'cleanning.avi';
workingDir = 'prototype\slides\results\detecting cells';
createVideo(workingDir , videoName);
toc