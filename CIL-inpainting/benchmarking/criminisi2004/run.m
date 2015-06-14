% Tabula rasa.
clear all
close all
clc

% Create mex file.
disp 'Compiling bestexemplarhelper...'
mex bestexemplarhelper.c

% Run example
disp 'Be patient! Calculating inpainting...'
[i1,i2,i3,c,d,mov] = inpaint_criminisi('data/bungee.png','data/bungee_mask.png');
disp 'Done!'
plotall;
pause
close all;
disp 'Showing evolution of inpainting...'
movie(mov);