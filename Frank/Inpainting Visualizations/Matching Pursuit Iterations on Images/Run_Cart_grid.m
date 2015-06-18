% __ This script is a visualization tool used in gaining better intuition
% of the matching pursuit algorithm. The intent is to create better than
% greedy algorithms for the purpose of inpainting.
close all; clear;

% Names of folders
ImageFolder = '\Image_Sources\';
DictionaryFolder = '\Dictionaries\';
HelperFunctionsFolder = '\Sub_functions\';
CurrentDir = pwd();
k = 1;

Errors = []; % mean squared errors for each image would be stored here
% Get all the Image related files
img_file_list = dir(strcat(CurrentDir,ImageFolder));
 
% Add extra functions and files to current directory
addpath(strcat(CurrentDir,ImageFolder))
addpath(strcat(CurrentDir,HelperFunctionsFolder))
addpath(strcat(CurrentDir,DictionaryFolder))
%% Parameters are inserted here
d = 16; %block size
sigma = 0.005; %stopping criterion for sparsecoding
rc_min = 0.03; %stopping criterion for sparsecoding
 
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
    
    mask=ones(size(I));
    mask = imread(mask_name);
    mask=mask./255;
    mask=mask(:,:,1);
    
%     mask=ones(size(I));
%     for i=1:size(mask,1)*size(mask,2)/1.01
%         a=floor(rand(1)*size(mask,1)+1);
%         b=floor(rand(1)*size(mask,2)+1);
%         mask(a,b)=0;
%     end
%     
    I_mask = I;
    I_mask(~mask) = 0;
    
    % Vectorize image matrix
    [I_vec, Im_size]= extract(I_mask,d);
%     [I_rec, ~]= extract(I,d);
%     [mask, ~]= extract(mask,d);
    
    % Get dictionary 
    
    U = buildDictionary(d^2,0);  %    0 second argument for CUSTOM Dictionary
    
%     [U1,Z, U] = dictionary_learning(I_rec,Im_size)
    
    % run modified visualization sparsecoding code
    t =cputime;
    I_out = C0nfidence_first( I_mask, mask); %This file here is a bit modified
    cputime - t
    
    Errors(k) = mean(mean(mean( ((I - I_out) ).^2)));
    
    k = k+1;
end


Result(1) = mean(Errors);

disp(['Average quadratic error: ' num2str(Result(1))])