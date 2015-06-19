close all
clear

dataDir = '../output/Dictionary efficiency';
outputDir = [dataDir,'/figures'];

fileList = dir([dataDir, '/*.mat']); 

for i=1:length(fileList)
    filename = fullfile(dataDir,fileList(i).name);
    load(filename);
    
    % Fill in the fields of the mat.
    allErrors(i) = {errors};
    allRuntimes(i) = {runtime};
    allSparsities(i) = {sparsity};

    % Name for the legend
    legendfullname = fileList(i).name;
    legendname(i) = {legendfullname(1:end-4)};
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
    figure(1);
    hold on
    scatter(allRuntimes{i}, allErrors{i}, 'MarkerFaceColor', colors(i),'MarkerEdgeColor','k');
    legend(legendname)
end

for i=1:length(allErrors)
    figure(1);
    linespec=['-','.',colors(i)];
    hold on
    errorbar(mean(allRuntimes{i}), mean(allErrors{i}), std(allErrors{i}) ,linespec, 'LineWidth',2);
    figure(1);
    hold on
    herrorbar(mean(allRuntimes{i}), mean(allErrors{i}), std(allRuntimes{i}) ,linespec);
end

grid on
xlabel('Runtime (s)');
yl = ylim(gca);
xl = xlim(gca);
xlim([0, xl(2)]);
%y = ylabel('Runtime (s)', 'fontsize', fontsize, 'rot', 0);
%set(y, 'position', [-0.22,yl(2)+0.0001]);
y = ylabel('Avg. quadratic Error');
axis([0 100 -0.002 0.01]);
title(dataDir(11:end))

%set(gca,'FontSize',fontsize)
%set(findall(gcf,'type','text'),'FontSize', fontsize)

saveas(gcf, fullfile(outputDir, 'scatter_plot'), 'pdf');