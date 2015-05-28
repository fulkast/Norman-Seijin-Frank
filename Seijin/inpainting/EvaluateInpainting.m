% Measure approximation error and compression ratio for several images.
%
% NOTE Images must be have .png ending and reside in the same folder.

file_list = dir(); 
k = 1;

Errors = []; % mean squared errors for each image would be stored here


for i = 3:length(dir) % running through the folder
    
    file_name = file_list(i).name; % get current filename
    
    % Only keep the images in the loop
    if (length(file_name) < 5)
        continue;
    elseif ( max(file_name(end-4:end) ~= '2.png'))
        continue;
    end
    mask_name = [file_name(1:end-4) '_mask.png'];
        
    % Read image, convert to double precision and map to [0,1] interval
    I = imread(file_name); 
    I = double(I) / 255; 
    I=I(1:512,1:512,1);
    % Read the respective binary mask
    % EVALUATION IS DONE WITH A FIXED MASK
    mask=ones(size(I));
    mask = imread(mask_name);
    mask=mask;
   
    
    %pepper noise
    %mask=ones(size(I));
    for i=1:size(mask,1)*size(mask,2)/2*0
        a=floor(rand(1)*size(mask,1)+1);
        b=floor(rand(1)*size(mask,2)+1);
        mask(a,b)=0;
    end

 mask=mask(:,:,1);
    I_mask = I;
    I_mask(~mask) = 0;
mask=im2double(mask);
imshow(I_mask)
    % Call the main inPainting function
   I_rec = inPainting(I_mask, mask);
   imshow(I_rec);
    % Measure approximation error
    Errors(k) = mean(mean(mean( ((I - I_rec) ).^2)));
    disp(Errors(k));
    k = k+1;
    break;
end

Result(1) = mean(Errors);

disp(['Average quadratic error: ' num2str(Result(1))])
