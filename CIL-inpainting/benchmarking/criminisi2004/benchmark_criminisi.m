function benchmark_criminisi(inputDir, outputDir, pathToMask)
% benchmark_criminisi runs in-painting on the images provided in input dir.
% The implementation is based on the proposition by Criminisi et. al 2004,
% and is provided by Bhat 
%
% INPUT
% inputDir:     Contains files according to the following pattern:
%               image       <filename>.png 
%               mask_image  <filename>_mask.png
%               image and mask_image must be of same size.
%
% outputDir:    Folder to which the results will be written.
%
% pathToMask:   Optional. If set, this mask will be used for all images
%               If not, the above file pattern will be applied.
%
% OUTPUT:       The script writes result into outputDir:
%               image_rec   <filename>_out.png

% Constants.
OUT_SUFFIX  = '_out';
MASK_SUFFIX = '_mask';

useSingleMask = false;
if nargin == 3
    useSingleMask = exist(pathToMask, 'file');
    maskPath = pathToMask;
end 

if ~exist(inputDir, 'dir')
    error(['Input directory ', inputDir, ' does not exist.'])
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

if ~exist(outputDir, 'dir')
    error(['Output directory ', outputDir, ' could not be created.']);
end

if isempty(dir('*.mex*'))
    disp 'Compiling bestexemplarhelper...'
    mex bestexemplarhelper.c
end

fileList = dir(inputDir); 


for i = 3:length(fileList) 
    
    % Input image and mask.
    imageName = fileList(i).name;
    [~, imageName, imageExt] = fileparts(imageName);
    imagePath = fullfile(inputDir, [imageName, imageExt]);
    if ~useSingleMask
        maskPath = fullfile(inputDir, [imageName MASK_SUFFIX imageExt]);
    end
    
    % Skip file if...
    %   - mask does not exist for this file
    %   - imageName contains "_mask"
    %   - it is not a png.
    if ~(exist(imagePath, 'file') && exist(maskPath, 'file')) ... 
            || ~isempty(strfind(imageName, '_mask')) ...
            || strcmp(imageExt, '.png' == 0)
        continue;
    end
    
    fprintf('Processing image %s...\n', imagePath);
    
    % Run inpaint script
    % TODO: speedup trick (because it takes ages otherwise)
    % TODO: check if it works with our type of mask
    % TODO: save output image
    % TODO: return statistics
    
    I_rec = inpaint_criminisi(imagePath, maskPath);
    I_rec = uint8(I_rec);
    
    % Save image to output dir.
    outPath = fullfile(outputDir, [imageName OUT_SUFFIX imageExt]);
    imwrite(I_rec, outPath);
end