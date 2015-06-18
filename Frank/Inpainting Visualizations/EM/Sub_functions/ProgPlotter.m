function ProgPlotter(atomIndex,U_thedictionary,...
    Z_all_coefficients,X_fullMatrix,InnerProduct_vector,PatchSideLength)


        subplot(1,3,1) 
imshow(imresize(reshape...
(U_thedictionary*Z_all_coefficients(:,atomIndex),PatchSideLength,PatchSideLength),16))
% Note the 16 above is just for magnifying the image and is not there
% because of the specific patch size
        title('Current Progress')
        subplot(1,3,2) 
imshow(imresize(reshape(X_fullMatrix(:,atomIndex),PatchSideLength,PatchSideLength),16))
        title('Actual Patch Trying to Fit')
        subplot(1,3,3) 
plot(sort(InnerProduct_vector))
        title('Coherence of Residue to Each Atom (sorted in ascending order)')
        xlabel('Atoms Sorted')
        ylabel('Coherence Value')
        
        drawnow;