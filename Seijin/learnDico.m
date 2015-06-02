
dataDir = 'dataPics';
file_list = dir(dataDir); 
k = 1;
siz=512;
Errors = []; % mean squared errors for each image would be stored here
IM=zeros(siz);

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
    if (ceil(siz/min(a,b))>1)
        I=expand(I,ceil(siz/min(a,b)));
    elseif (floor(min(a,b)/siz)>1)
        I=compress(I,floor(min(a,b)/siz));
    end
    
    I=I(1:siz,1:siz,1);
    IM=[IM I];
end

X=my_im2col(IM,16,0,[0 0]);
dictionary_learning(X,size(IM));