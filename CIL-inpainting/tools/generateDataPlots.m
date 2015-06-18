dataDir = '../output';
outputDir = '../output/figures';

fileList = dir([dataDir, '/*.mat']); 

for i=1:length(fileList)
    filename = fullfile(dataDir,fileList(i).name);
    load(filename);
    
    % Fill in the fields of the mat.
    allErrors(i) = {errors};
    allRuntimes(i) = {runtime};
    allSparsities(i) = {sparsity};
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Scatter plot
close all;

fontsize = 15;
colors = ['b', 'g', 'r', 'c', 'm', 'k'];
%set(0,'DefaultAxesFontName', 'Gill Sans MT')
set(0,'DefaultAxesFontName', 'Times New Roman')
%set(0,'DefaultAxesFontName', 'Helevtica')
set(0,'DefaultAxesFontSize', 14)

for i=1:length(allErrors)
    figure(1); scatter(allRuntimes{i}, allErrors{i}, 'MarkerFaceColor', colors(i));
end

grid on
xlabel('Avg. quadratic Error');
yl = ylim(gca);
xl = xlim(gca);
xlim([0, xl(2)]);
%y = ylabel('Runtime (s)', 'fontsize', fontsize, 'rot', 0);
%set(y, 'position', [-0.22,yl(2)+0.0001]);
y = ylabel('Runtime (s)');
title('Speed ')

%set(gca,'FontSize',fontsize)
%set(findall(gcf,'type','text'),'FontSize', fontsize)

saveas(gcf, fullfile(outputDir, 'scatter_plot'), 'pdf');