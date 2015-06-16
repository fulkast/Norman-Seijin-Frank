
dataDir = 'DataPics';
file_list = dir(dataDir); 
k = 1;
siz=512;
Errors = []; % mean squared errors for each image would be stored here

N=0;
for i = 3:length(file_list) % running through the folder
    
    file_name = fullfile(dataDir, file_list(i).name); % get current filename
    
   % file_name = file_list(i).name; % get current filename
    
    % Only keep the images in the loop
    if ( (file_name(end-7:end) == 'mask.png'))
        continue;
    end
    try
    I = imread(file_name); 
   catch
       continue;
   end
    N=N+1;
end

IM=zeros(siz*sqrt(N));
Nn=0;
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
    [a b]=size(I);
    if (ceil(siz/min(a,b))>1)
        I=expand(I,ceil(siz/min(a,b)));
    elseif (floor(min(a,b)/siz)>1)
        I=compress(I,floor(min(a,b)/siz));
    end
    
    I=I(1:siz,1:siz,1);
    IM(mod(Nn,sqrt(N))*siz+1:mod(Nn,sqrt(N))*siz+siz,floor(Nn/sqrt(N))*siz+1:floor(Nn/sqrt(N))*siz+siz)=I(:,:);
    Nn=Nn+1;
end

X=my_im2col(IM,16,0,[0 0]);
% save('DATA.mat','X');
% Y=X;
% Y(1,:)=0.*Y(1,:);
% imshow((my_col2im(Y,16,size(IM),0,[0 0]))  );
imshow(IM);
figure
dictionary_learning(X,size(IM));