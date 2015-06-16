%Fixing Neighbour Plotter

function NeibourIntrudingPlotter(InitialPatch,Uloc,...
    Zloc,PatchSideLength)


        subplot(1,3,1) 
imshow(imresize(reshape...
(Uloc*Zloc,PatchSideLength,PatchSideLength),16))
% Note the 16 above is just for magnifying the image and is not there
% because of the specific patch size
        title('Current Progress')
        subplot(1,2,2) 
imshow(imresize(reshape(InitialPatch,PatchSideLength,PatchSideLength),16))
        title('Actual Patch Trying to Fit')
        
        drawnow;