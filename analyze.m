initCobraToolbox
changeSolver('gurobi')
load('recon1.mat');
modelGenes = metabolicmodel.genes;

% Get gene symbols that are differentially expressed based on proteomics
% data
%geneExp = readcell('DestackProteomics2010.xlsx');
[num, txt] = xlsread('DestackProteomics2010.xlsx');
txt(1, :) = [];

upPos = num(:, 7) > 2;
downPos = num(:, 7) < -2;
upGenes = txt(upPos);
downGenes = txt(downPos);
upGenes = intersect(upGenes, modelGenes);
downGenes = intersect(downGenes, modelGenes);


