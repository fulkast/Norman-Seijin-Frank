
% The training Data directory. It should contain a square number of
% 512x512 pictures
dataDir = 'TrainingData';
file_list = dir(dataDir); 
siz=512;
N=0;
for i = 3:length(file_list) % running through the folder to count the number of images
    file_name = fullfile(dataDir, file_list(i).name); 
    try
    I = imread(file_name); 
   catch
       continue;
   end
    N=N+1;
end
% IM is the concatenation image of all training pictures
IM=zeros(siz*sqrt(N));
Nn=0;
for i = 3:length(file_list) % running through the folder
    file_name = fullfile(dataDir, file_list(i).name); % get current filename
    fprintf('Reading file %s...\n', file_name)
   try
    I = imread(file_name); 
   catch
       continue;
   end
    I = double(I) / 255; 
    %Concatenating images
    IM(mod(Nn,sqrt(N))*siz+1:mod(Nn,sqrt(N))*siz+siz,floor(Nn/sqrt(N))*siz+1:floor(Nn/sqrt(N))*siz+siz)=I(:,:);
    Nn=Nn+1;
end

X=my_im2col(IM,16,0,[0 0]);
imshow(IM);
figure
dictionary_learning(X,size(IM));