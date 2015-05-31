% Tabula rasa.
clear all
close all
clc

% Create mex file.
disp 'Compiling bestexemplarhelper...'
mex bestexemplarhelper.c

% Run example
disp 'Be patient! Calculating inpainting...'
[i1,i2,i3,c,d,mov] = inpaint('bungee0.png','bungee1.png',[0 255 0]);
disp 'Done!'
plotall;
pause
close all;
disp 'Showing evolution of inpainting...'
movie(mov);