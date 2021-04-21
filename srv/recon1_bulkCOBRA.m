%% RECON1 A549 Bulk COBRA Analysis
% *Author*: Scott Campit
%% Summary
% This script computes the reaction knockouts and metabolic fluxes for the scCOBRA 
% models.
%% 1. Load up COBRA Toolbox
% This code block loads up the COBRA Toolbox for constraint-based metabolic 
% modeling.

% Initialize metabolic modeling components
clear all;

% Needed to load COBRA onto Great lakes
addpath('/nfs/turbo/umms-csriram/scampit/Software/cobratoolbox');

initCobraToolbox(false); changeCobraSolver('gurobi', 'all');
%% 2. Set flux balance analysis hyperparameters
% Assign default values of parameters for the iMAT algorithm.

kap = [];       rho = [];
isgenes = true; pfba = true;
eps = [];       eps2 = [];
%% 3. Compute COBRA data
% This code block fits the metabolic reconstruction with differentially expressed 
% genes using the iMAT algorithm (|constrain_flux_regulation.m|).
% A. Setting up parallel computations
% We can set the number of CPUs to use for our calculations.
% 
% The code below is for running this script on Great Lakes

% Initialize the Umich Cluster profiles
%setupUmichClusters

% We get from the environment the number of processors
%NP = str2num(getenv('SLURM_NTASKS'));

% Create the pool for parfor to use
%thePool = parpool('current', NP);
%poolobj = parpool;
%addAttachedFiles(poolobj, {})
%% 
% This codeblock is for running this script locally.

%workers = 4;
%parpool("local", workers);
parpool
% B. Perform bulk simulations
% This code block will perform flux balance analysis and knockout analysis for 
% the three bulk datasets.
% i. Garcia
% This is the EMT proteomics dataset from the Garcia lab.

% DELL
%load D:/Analysis/EMT/data/Garcia.mat

% NFS
load           ~/Turbo/scampit/Analysis/EMT/data/Garcia.mat
garcia_path = "~/Turbo/scampit/Analysis/EMT/garcia/recon1/";
%% 
% First, extract the datasets.

entrez_zdata = data;
entrez       = geneids;
%% 
% Then run the flux and knockout simulations.

filenames = ["fc05min.mat", "fc60min.mat", ...
             "fc24hr.mat", "fc48hr.mat"];

recon1.genes = recon1.geneEntrezID;

parfor j = 1:size(entrez_zdata, 2)
    % Grab indicies for significant differentially regulated genes
    cell_data = entrez_zdata(:, j);
    up_idx    = cell_data > 0;
    down_idx  = cell_data < 0;
        
    % Get up- and down-regulated gene identifiers
    upgenes    = entrez(up_idx);
    downgenes  = entrez(down_idx);
    
    % Force rho to be 10
    rho = repelem(10, length(upgenes));
    
    % Fit models and get data
    [soln, ~, ~, cell_mdl] = constrain_flux_regulation(recon1, ...
                                                       upgenes, downgenes, ...
                                                       kap, rho, ...
                                                       eps, ...
                                                       isgenes, ...
                                                       pfba);
    [geneKO, rxnKO]  = knockOut(cell_mdl, 'All');
    
    % Save data as individual files
    cobra_cell        = struct();
    cobra_cell.id     = Garcia.fcid{j};
    cobra_cell.geneko = geneKO;
    cobra_cell.rxnko  = rxnKO;
    cobra_cell.flux   = soln;
    parsave(strcat(garcia_path, filenames(j)), cobra_cell, j);   
end
% ii. GSE17518
% This is EMT RNASeq data from Thannickal et al.

% DELL
%load D:/Analysis/EMT/data/GSE17518.mat

% NFS
load             ~/Turbo/scampit/Analysis/EMT/data/GSE17518.mat
gse17518_path = "~/Turbo/scampit/Analysis/EMT/gse17518/recon1/";
%% 
% First, extract the datasets

entrez_zdata = data;
entrez       = geneids;
%% 
% Then run the flux and knockout simulations.

cell_data = entrez_zdata;
up_idx    = cell_data > 0;
down_idx  = cell_data < 0;
    
% Get up- and down-regulated genes
upgenes    = entrez(up_idx);
downgenes  = entrez(down_idx);

% Force rho to be 10
rho = repelem(10, length(upgenes));

% Fit models and get 
[soln, ~, ~, cell_mdl] = constrain_flux_regulation(recon1, ...
                                                   upgenes, downgenes, ...
                                                   kap, rho, ...
                                                   eps, ...
                                                   isgenes, ...
                                                   pfba);
[geneKO, rxnKO]  = knockOut(cell_mdl, 'All');

% Save data as individual files
cobra_cell        = struct();
cobra_cell.id     = GSE17518.fcid{1};
cobra_cell.geneko = geneKO;
cobra_cell.rxnko  = rxnKO;
cobra_cell.flux   = soln;

j = 1;
parsave(strcat(gse17518_path, "fc72.mat"), cobra_cell, j);   
% iii. GSE17708
% Now let's perform some knockouts. I will perform gene knockouts to save time.

% DELL
%load D:/Analysis/EMT/data/GSE17708.mat

% NFS
load             ~/Turbo/scampit/Analysis/EMT/data/GSE17708.mat
gse17708_path = "~/Turbo/scampit/Analysis/EMT/gse17708/recon1/";
%% 
% First, extract the datasets

entrez_zdata = data;
entrez       = geneids;
%% 
% Then run the flux and knockout simulations.

filenames = ["fc05min.mat", "fc01hr.mat", "fc02hr.mat", ...
             "fc04hr.mat", "fc08hr.mat", "fc16hr.mat", ...
             "fc24hr.mat", "fc72hr.mat"];

parfor j = 1:size(entrez_zdata, 2)
    cell_data = entrez_zdata(:, j);
    up_idx    = cell_data > 0;
    down_idx  = cell_data < 0;
        
    % Get up- and down-regulated genes
    upgenes    = entrez(up_idx);
    downgenes  = entrez(down_idx);
    
    % Force rho to be 10
    rho = repelem(10, length(upgenes));
    
    % Fit models and get 
    [soln, ~, ~, cell_mdl] = constrain_flux_regulation(recon1, ...
                                                       upgenes, downgenes, ...
                                                       kap, rho, ...
                                                       eps, ...
                                                       isgenes, ...
                                                       pfba);
    [geneKO, rxnKO]  = knockOut(cell_mdl, 'All');
    
    % Save data as individual files
    cobra_cell        = struct();
    cobra_cell.id     = GSE17708.fcid{j};
    cobra_cell.geneko = geneKO;
    cobra_cell.rxnko  = rxnKO;
    cobra_cell.flux   = soln;
    parsave(strcat(gse17708_path, filenames(j)), cobra_cell, j);   
end
%% Summary
% This notebook goes through bulk and single cell flux balance analysis. To 
% get differentially active or differentially sensitive metabolic reactions, you 
% can check out the |scDRSA.mlx| file.