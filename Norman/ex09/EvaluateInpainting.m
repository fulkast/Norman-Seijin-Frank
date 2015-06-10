% Measure approximation error and compression ratio for several images.
%
% NOTE Images must be have .png ending and reside in the same folder.
dataDir = 'data';
file_list = dir(dataDir); 
k = 1;

Errors = []; % mean squared errors for each image would be stored here


for i = 3:length(file_list) % running through the folder
    
    file_name = fullfile(dataDir, file_list(i).name); % get current filename
    
    % Only keep the images in the loop
    if (length(file_name) < 5)
        continue;
    elseif ( max(file_name(end-4:end) ~= '2.png'))
        continue;
    end
    mask_name = [file_name(1:end-4) '_mask.png'];
        
    % Read image, convert to double precision and map to [0,1] interval
    fprintf('Reading file %s...\n', file_name)
    I = imread(file_name); 
    I = double(I) / 255; 
    
    % Read the respective binary mask
    % EVALUATION IS DONE WITH A FIXED MASK
    mask = imread(mask_name);
    %mask = random_mask(size(I,1), 0.7);
    
    I_mask = I;
    I_mask(~mask) = 0;
          
    % Call the main inPainting function
    fprintf('Running the in-pating algorithm...\n');
    
    sf = 1.;
    I = imresize(I, sf, 'bilinear');
    I_mask = imresize(I_mask, sf, 'nearest');
    mask = imresize(mask, sf, 'nearest');
    
    I_rec = inPainting(I_mask, mask);
    
    figure(1)
    imshow(I_mask);
    figure(2)
    imshow(I_rec);
    pause;
    
    % Measure approximation error
    Errors(k) = mean(mean(mean( ((I - I_rec) ).^2)));
    
    k = k+1;
end

Result(1) = mean(Errors);

disp(['Average quadratic error: ' num2str(Result(1))])
