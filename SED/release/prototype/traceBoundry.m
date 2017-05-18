workingDir = 'prototype\slides\results\flow';

imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
a = size(imageNames);
length = a(1);

slide = imread(fullfile(workingDir,'382.png'));

workingDir = 'prototype\slides\results\SED + Sobel';
nextSlide = imread(fullfile(workingDir,'383.png'));

imshow(nextSlide)

[B,L] = bwboundaries(slide,'noholes');
%imshow(label2rgb(L, @jet, [.5 .5 .5]))

boundary = B{1};
newBoundary = boundary;


%%

hold on
plot(boundary(:,2), boundary(:,1),...
   'r','LineWidth',2);

newPoly = extendPoly(boundary,15);
newBoundary2 = newPoly{1,1};

hold on
plot(int16(newBoundary2(:,2)), int16(newBoundary2(:,1)),...
   'g','LineWidth',2);
%%
mat = zeros(size(slide));
for ii=1:1105
    mat(int16(newBoundary2(ii,1)), int16(newBoundary2(ii,2)))= 1; 
end

[i,j]=find(mat); % find i,j coordinates of non-zero elements
k=convhull(i,j); % create convexhull
[iv,jv]=ndgrid(1:size(mat,1),1:size(mat,2)); % create index space
mat = inpolygon(iv,jv,i(k),j(k)) > 0; % create filled-in matrix:
imshow(mat)

%%
nextSlide(mat == 0 ) = 255;
imshow(nextSlide)

%% Active contur


bw = activecontour(nextSlide,mat,300, 'edge');

figure
imshow(bw)
title('Segmented Image')


%%
slide = nextSlide;

workingDir = 'prototype\slides\results\SED + Sobel';
nextSlide = imread(fullfile(workingDir,'384.png'));

imshow(nextSlide)

[B,L] = bwboundaries(slide,'noholes');
%imshow(label2rgb(L, @jet, [.5 .5 .5]))

boundary = B{1};
newBoundary = boundary;
%%

hold on
plot(boundary(:,2), boundary(:,1),...
   'r','LineWidth',2);

newPoly = extendPoly(boundary,15);
newBoundary2 = newPoly{1,1};

hold on
plot(int16(newBoundary2(:,2)), int16(newBoundary2(:,1)),...
   'g','LineWidth',2);

%%
workingDir = 'prototype\slides\results\flow';
slide = imread(fullfile(workingDir,'382.png'));
stratFrame = 382;

workingDir = 'prototype\slides\results\SED + Sobel';
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
a = size(imageNames);
len = a(1);

for ii=stratFrame+1:408
    nextSlide = imread(fullfile(workingDir,[num2str(ii), '.png']));
    nextSlide = nextSlide < 5;
    
    grenn_slice = imdilate(slide ,strel('disk' , 10));
    grenn_slice = imfill(grenn_slice, 'holes');
    %nextSlide = 255 - nextSlide;
    nextSlide(grenn_slice == 0 ) = 0;
    slide = nextSlide;
    imshow(nextSlide)
    pause(0.2);
end

