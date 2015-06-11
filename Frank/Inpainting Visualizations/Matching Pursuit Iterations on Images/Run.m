% __ This script is a visualization tool used in gaining better intuition
% of the matching pursuit algorithm. The intent is to create better than
% greedy algorithms for the purpose of inpainting.
close all; clear;

% Names of folders
ImageFolder = '\Image_Sources\';
DictionaryFolder = '\Dictionaries\';
HelperFunctionsFolder = '\Sub_functions\';
CurrentDir = pwd();

% Get all the Image related files
img_file_list = dir(strcat(CurrentDir,ImageFolder));

% Add extra functions and files to current directory
addpath(strcat(CurrentDir,ImageFolder))
addpath(strcat(CurrentDir,HelperFunctionsFolder))
addpath(strcat(CurrentDir,DictionaryFolder))
%% Parameters are inserted here
d = 16; %block size
sigma = 0.01; %stopping criterion for sparsecoding
rc_min = 0.01; %stopping criterion for sparsecoding
 
%%
for i = 3:length(img_file_list) % runing through the folder
    
    file_name = img_file_list(i).name;

    % Only keep the images in the loop
    if (length(file_name) < 5)
        continue;
    elseif ( max(file_name(end-4:end) ~= '2.png')) % this bit assumes that the image ends in 512.png
                                                    % and the mask ends in
                                                    % something else (shown
                                                    % below)
        continue;
    end
    mask_name = [file_name(1:end-4) '_mask.png'];
        
    % Read image, convert to double precision and map to [0,1] interval
    I = imread(file_name); 
    I = double(I) / 255;     
    
    % Vectorize image matrix
    [I_vec, Im_size]= extract(I,d);
    
    % Get dictionary
    
    U = buildDictionary(d^2); 
    
    % run modified visualization sparsecoding code
    [Z] = sparseCoding4visual(U, I_vec, ones(size(I_vec)), sigma, rc_min, Im_size); %This file here is a bit modified
    
    
    
    
    
    
end