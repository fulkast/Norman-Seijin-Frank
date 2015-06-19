%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main execution script:
% Please check README before running
%

% Path to data dir.
dataDir = 'data/selection';

% Type of mask:
%   0: reads the mask from a file
%   1: create salt and pepper noise
maskType = 0;

% The amount of noise added to the image.
% Must be a value between 0 and 1.
% Has an effect only for maskType == 0
maskPepperness = 0.7;

% Show images and pause execution.
showResults = false;

% Use fixed mask image (e.g. mask.png). This option is useful if not every
% input images comes with its own mask.
useFixedMaskFile = true;

% Path to the fixed mask image. Use only if useFixedMaskFile == true
fixedMaskFilePath = fullfile(dataDir, 'mask.png');

% Set useSpecialDictionary if a special dictionary should be copied into
% the working directory. The suffix will be used when the runtime
% statistics are saved to a .mat file.
%
% !!!WARNING!!! If useSpecialDictionary==true, the existing
%               dictionary.mat will be overwritten!!!
useSpecialDictionary = true;
dictionaryFilePath = 'dictionaries/projectedKmeanHS16pics.mat';  
dictionarySuffix = 'projectedKmeanHS16pics';

% Set true to save the resulting images and statistics to outputDir.
saveResults = true;

% Path to output directory. The folder will be created if it doesn't exist
% yet. Existing content might be overwritten.
outDir = 'output';


%% Setup
fileList = dir(dataDir); 

% Mean squared errors.
errors = []; 
count = 1;

% Control random number generator (in case it is used).
rng('default');
rng(1);

% Copy dictionary to dictionary.mat
if useSpecialDictionary
    targetDictionary = 'dictionary.mat';
    if exist(dictionaryFilePath, 'file') && ~isequal(dictionaryFilePath, targetDictionary)
        fprintf('Copying dictionary %s to %s\n', dictionaryFilePath, targetDictionary);
        copyfile(dictionaryFilePath, targetDictionary, 'f');
    else
        error('Something went wrong. Check dictionaryFilePath.')
    end
end

%% Loop over data dir
for i = 3:length(fileList) 
    
    % Input image and mask.
    fileName = fileList(i).name;
    filePath = fullfile(dataDir, fileName);
    if (strcmp(filePath(end-7:end), 'mask.png'))
        continue;
    end
    if useFixedMaskFile
        maskName = fixedMaskFilePath;
    else
        maskName = [filePath(1:end-4) '_mask.png'];
    end
        
    % Read image.
    fprintf('Reading file %s...\n', filePath);
    try
        I = imread(filePath);
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
            mask = imread(maskName);
        case 1
            mask = random_mask(size(I), maskPepperness);
    end
    mask=im2double(mask);
    
    % Apply mask.
    I_mask = I;
    I_mask(~mask) = 0;
    
    % In-painting.
    [I_rec, stats] = inPainting(I_mask, mask);
    runtime(count) = stats.runtime;
    sparsity(count) = stats.sparsity;
    
    if saveResults
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        outFile = fullfile(outDir, fileName);
        fprintf('Saving file %s...\n', outFile); 
        imwrite(I_rec, outFile);
    end
    
    % Measure approximation error and print result.
    errors(count) = mean(mean(mean(((I - I_rec) ).^2)));
    fprintf('error= %s  \n',errors(count));
    count = count + 1;
    
    if showResults
        figure(1)
        imshow(I_mask);
        figure(2)
        imshow(I_rec);
        pause
    end
end

%%
result(1) = mean(errors);
disp(['Average quadratic error: ' num2str(result(1))])

if saveResults
    save(fullfile(outDir, ['stats_', dictionarySuffix, '.mat']), 'errors', 'runtime', 'sparsity');
end
