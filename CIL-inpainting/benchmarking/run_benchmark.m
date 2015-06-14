% Execute benchmark.
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
% verbose:      Print additional information to console.
%
% OUTPUT:       The script writes result into outputDir:
%               image_rec   <filename>_out.png

%% Settings .
inputDir = 'data';
outputDir = 'output';
pathToMask = '';
verbose = true;
runCriminisi = false;
runBIG = true;
runRef = true;

% Mask type
%   0: mask from file
%   1: salt'n'pepper noise
maskType = 0;

% Path to single mask (if not every image comes with it's own mask)
pathToMask = '';

useSingleMask = exist(pathToMask, 'file');
maskPath = pathToMask;
    
%% Constants.
PATH_CRIMINISI = 'criminisi2004/';
PATH_BIG = '..';

addpath(PATH_CRIMINISI);
addpath(PATH_BIG);

% Image suffixes.
IN_SUFFIX   = '_in';
OUT_SUFFIX  = '_out';
MASK_SUFFIX = '_mask';

% Method suffixes.
CRIM_SUFFIX = '_crim';
BIG_SUFFIX = '_big';
REF_SUFFIX = '_ref';

% Mask suffix.
switch maskType
    case 0
        MASK_SUFFIX = '_mask';
    case 1
        MASK_SUFFIX = '_pepper';
end

%% Check inputs.
if ~exist(inputDir, 'dir')
    error(['Input directory ', inputDir, ' does not exist.'])
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

if ~exist(outputDir, 'dir')
    error(['Output directory ', outputDir, ' could not be created.']);
end

%% Check Criminisi
if isempty(dir([PATH_CRIMINISI, '*.mex*']))
    disp 'Compiling bestexemplarhelper...'
    mex('-outdir', PATH_CRIMINISI, fullfile(PATH_CRIMINISI, 'bestexemplarhelper.c'));
end

%% Run
count = 1;
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
    info.file = imagePath;
    info.mask = maskType;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Criminisi method.
    if runCriminisi
        starttime = cputime;
        [I_rec, ~, I_mask] = inpaint_criminisi(imagePath, maskPath, verbose);
        info.runtimeCriminisi = cputime - starttime;
        I_rec = uint8(I_rec);

        % Save images to output dir.
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX CRIM_SUFFIX imageExt]);
        imwrite(I_rec, outPath);
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX IN_SUFFIX imageExt]);
        imwrite(I_mask, outPath);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run BIG method.
    
    %TODO: unify interfaces, read images only once.
    
    % Read images.
    try
        I = imread(imagePath);
        if size(I,3) > 1
            I = rgb2gray(I);
        end
    catch
        continue;
    end
    
    % Convert to double precision and map to [0,1] interval.
    I = double(I) / 255;
    
    % Read mask.
    switch maskType
        case 0
            mask = imread(maskPath);
        case 1
            mask = random_mask(size(I), maskPepperness);
    end
    mask=im2double(mask);
    
    % Apply mask.
    I_mask = I;
    I_mask(~mask) = 0;
    
    % Settings
    settings.overlap=8;         % even, less or equal to half the size of neib
    settings.rc_min = 0.01;     % rc_min: minimal residual correlation before stopping
    settings.neib = 16;         % neib: The patch sizes used in the decomposition of the image
    settings.sigma = 0.01;      % sigma: residual error stopping criterion, normalized by signal norm
    settings.verbose = verbose; %
    
    if runBIG
        starttime = cputime;
        I_rec = inPainting(I_mask, mask, settings);
        info.runtimeBIG = cputime - starttime;

        % Save images to output dir.
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX BIG_SUFFIX imageExt]);
        imwrite(I_rec, outPath);
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX IN_SUFFIX imageExt]);
        imwrite(I_mask, outPath);
    end
    
    if runRef
        % Disable overlapping patches.
        settings.overlap=0;
        
        starttime = cputime;
        I_rec = inPainting(I_mask, mask, settings);
        info.runtimeRef = cputime - starttime;

        % Save images to output dir.
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX REF_SUFFIX imageExt]);
        imwrite(I_rec, outPath);
        outPath = fullfile(outputDir, [imageName MASK_SUFFIX IN_SUFFIX imageExt]);
        imwrite(I_mask, outPath);
    end
    
    infos{count} = info;
    
    count = count + 1;
end