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


