%% Virabiales
rootDir = 'prototype\slides';
items = ['10000108_1_CTce_ThAb.nii.gz';'10000109_1_CTce_ThAb.nii.gz';'10000111_1_CTce_ThAb.nii.gz';'10000113_1_CTce_ThAb.nii.gz';'10000127_1_CTce_ThAb.nii.gz';];
    items = cellstr(items);
    volume = cell(length(items),1);
    sedSeg = cell(length(items),1);
    sobelOnSEDseg = cell(length(items),1);
    GT_slides = cell(length(items),1);
    finalSeg = cell(length(items),1);
    startFrame = cell(length(items),1);
    endFrame = cell(length(items),1);
    GTframe = [274,238,226,255,271];
%GTframe = [382,548,521,536,520];
%GTframe = [337,493,458,467,445];

%% Loading the CT scans.
tic
disp('Loading scans...');
for ii=1:length(items)
    CTscan = load_untouch_nii_gzip(strjoin(items(ii)));
    volume{ii} = CTscan.img;
    createCTslides(cell2mat(volume(ii)), strjoin(items(ii)), rootDir);
end
clearvars CTscan
toc

%% creating the SED segmentation of CT scan.
tic
disp('evaluating SED segmantation on scans...');
for ii=1:length(items)
    itemFolderName = strsplit(strjoin(items(ii)),'.');
    workingDir = fullfile(rootDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
    sedSeg{ii} = createSEDsegmentation(cell2mat(volume(ii)), workingDir);
end
clearvars volume
toc

%% creating the Sobel segmentation of SED segmentation.
tic
disp('evaluating Sobel segmantation on SED segmantations...');
for ii=1:length(items)
    itemFolderName = strsplit(strjoin(items(ii)),'.');
    workingDir = fullfile(rootDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
    sobelOnSEDseg{ii} = sobelOnSED(cell2mat(sedSeg(ii)), workingDir);
end
clearvars sedSeg
toc

%% get 1 slide ground truth segmentation
tic
disp('getting 1 slide ground truth segmentation...');
for ii=1:length(items)
    itemFolderName = strsplit(strjoin(items(ii)),'.');
    workingDir = fullfile(rootDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
   [ GT_slides{ii},startFrame{ii},endFrame{ii}] = getOneSlideGroundTruth(strjoin(items(ii)), GTframe(ii) , workingDir);
end
toc

%% trace the Liver for final segmentation
tic
disp('starting final segmentation...');
for ii=1:length(items)
    itemFolderName = strsplit(strjoin(items(ii)),'.');
    workingDir = fullfile(rootDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
    finalSeg{ii} = traceLiver( cell2mat(sobelOnSEDseg(ii)) ,cell2mat(GT_slides(ii)), cell2mat(startFrame(ii)), cell2mat(endFrame(ii)), GTframe(ii) ,workingDir );
end
clearvars sobelOnSEDseg
toc

%% saving the result
tic
disp('saving segmentation...');
for ii=1:length(items)
    itemFolderName = strsplit(strjoin(items(ii)),'.');
    workingDir = fullfile(rootDir,itemFolderName(1));
    workingDir = strjoin(workingDir);
    filename = strcat(strjoin(itemFolderName(1)) , '_58_result.nii.gz');
    fullname = fullfile(workingDir, filename);

    CTscan = load_untouch_nii_gzip(strjoin(items(ii)));
    CTscan.img = cell2mat(finalSeg(ii));
    save_untouch_nii_gzip(CTscan, fullname);
end
toc

