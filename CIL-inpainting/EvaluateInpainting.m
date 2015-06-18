%%
% Input
%
% maskType:
%   0: from file
%   1: salt and pepper noise
% maskPepperness
%   value between 0.0 or 1.0

dataDir = 'data';
maskType = 0;
maskPepperness = 0.7;
showResults = true;

outDir = 'output';
saveResults = false;

%% Setup
fileList = dir(dataDir); 

% Mean squared errors.
errors = []; 
count = 1;
siz=512;
% Control random number generator (in case it is used).
rng('default');
rng(1);
tic
%% Loop over data dir
for i = 3:length(fileList) 
    
    % Input image and mask.
    fileName = fileList(i).name;
    filePath = fullfile(dataDir, fileName);

    maskName = ['mask.png'];
        
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
    [a,b,c]=size(I);
    if (ceil(siz/min(a,b))>1)
        continue
        I=expand(I,ceil(siz/min(a,b)));
    elseif (floor(min(a,b)/siz)>1)
        I=compress(I,floor(min(a,b)/siz));
    end
    n1=floor(size(I,1)/2);
    n2=floor(size(I,2)/2);
    I=I(n1-siz/2+1:n1+siz/2,n2-siz/2+1:n2+siz/2,1);
    
    % Read mask.
    switch maskType
        case 0
            mask = imread(maskName);
        case 1
            mask = random_mask(size(I), maskPepperness);
    end
    mask=im2double(mask);
    I=I(1:512,1:512,1);
    % Apply mask.
    I_mask = I;
    I_mask(~mask) = 0;
    
    % In-painting.
    I_rec = inPainting(I_mask, mask);
    
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
%         pause
    end
end

result(1) = mean(errors);
disp(['Average quadratic error: ' num2str(result(1))])
toc
