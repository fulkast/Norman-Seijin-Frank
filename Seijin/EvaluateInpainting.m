% Measure approximation error and compression ratio for several images.
%
% NOTE Images must be have .png ending and reside in the same folder.
dataDir = 'Test set';
file_list = dir(dataDir); 
k = 1;
siz=512;
Errors = []; % mean squared errors for each image would be stored here

tic
for i = 3:length(file_list) % running through the folder
    
    file_name = fullfile(dataDir, file_list(i).name); % get current filename
    
   % file_name = file_list(i).name; % get current filename
    
    % Only keep the images in the loop
    if ( (file_name(end-7:end) == 'mask.png'))
        continue;
    end
   % mask_name = [file_name(1:end-4) '_mask.png'];
     mask_name = ['mask.png'];
        
    % Read image, convert to double precision and map to [0,1] interval
    fprintf('Reading file %s...\n', file_name)
   try
    I = imread(file_name); 
   catch
       continue;
   end
    I = double(I) / 255; 
    [a,b,c]=size(I);
    if (ceil(siz/min(a,b))>1)
        I=expand(I,ceil(siz/min(a,b)));
    elseif (floor(min(a,b)/siz)>1)
        I=compress(I,floor(min(a,b)/siz));
    end
    n1=floor(size(I,1)/2);
    n2=floor(size(I,2)/2);
    I=I(n1-siz/2+1:n1+siz/2,n2-siz/2+1:n2+siz/2,1);
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
%imshow(I_mask)
    % Call the main inPainting function
    
    tic
   I_rec = inPainting(I_mask, mask);
   toc
   %   figure
   imshow(I_rec);

    % Measure approximation error
    Errors(k) = mean(mean(mean( ((I - I_rec) ).^2)));
    fprintf('error= %s  \n',Errors(k));
    k = k+1;
   % break;
end
toc
Result(1) = mean(Errors);

disp(['Average quadratic error: ' num2str(Result(1))])
