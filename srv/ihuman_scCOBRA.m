%% Single Cell COBRA Analysis
% *Author*: Scott Campit
%% Summary
% This script computes the reaction knockouts and metabolic fluxes for the scCOBRA 
% models.
%% 1. Load Data
% First, let's load the COBRA toolbox.

% Initialize metabolic modeling components
clear all;
initCobraToolbox(false); changeCobraSolver('gurobi', 'all');
% Load iHUMAN
% This loads the iHUMAN metabolic reconstruction.

load D:/Data/Mappings/Reconstructions/Human1/ihuman_model_data_dec2020.mat ihuman_hams
model = ihuman_hams;
%% 
% Set the objective function to maximize biomass

% Set up objective function
model.c = zeros(size(model.c));
model.c(string(model.rxns) == 'biomass_human') = 1;
%% 2. Load gene expression dataset
% First, we'll load the A549 MAGIC-imputed gene expression dataset. 
% 
% For more information, see the |save_diffexp_data.mlx| livescript to see how 
% the initial dataset was preprocessed.

% DELL
%load D:/Analysis/EMT/scRNASeq_ihuman_diffexp.mat

% ACLX
%load ~/Analysis/EMT/scRNASeq_ihuman_diffexp.mat

% NFS
load ~/Turbo/scampit/Analysis/EMT/scRNASeq_ihuman_diffexp.mat
%% 3. Set flux balance analysis hyperparameters
% Assign default values of parameters for the iMAT algorithm.

eps = [];
kap = [];
rho = [];
isgenes = true;
eps2 = [];
pfba = true;
%% 4. Compute COBRA data
% This code block fits the metabolic reconstruction with differentially expressed 
% genes using the iMAT algorithm (|constrain_flux_regulation.m|).
% Setting up parallel computations
% We can set the number of CPUs to use for our calculations.

%workers = 4;
parpool("local", workers);
% iHUMAN metabolic reconstruction
% We'll use the iHUMAN metabolic reconstruction to perform flux balance analysis 
% and knockouts.
% 
% This is the save file path.

% ACLX
%filename = "~/Analysis/EMT/scRNASeq_ihuman_profiles.mat";

% NFS
filename = "~/Turbo/scampit/Analysis/EMT/scRNASeq_ihuman_profiles.mat";

% DELL
%filename = "D:/Analysis/EMT/scRNASeq_ihuman_profiles.mat";

% Save it
save(filename, 'scRxnKO', 'scRxnFx', 'scGeneKO');
% Single-cell simulations
% Now let's perform some knockouts. I will perform gene knockouts to save time.

scRxnFx  = zeros(length(time_course), length(model.rxns));
scGeneKO = zeros(length(time_course), length(model.genes));
scRxnKO  = zeros(length(time_course), length(model.rxns));

parfor j = 1:size(zdata, 2)
    
    cell_data = zdata(:, j);
    up_idx    = cell_data > 0;
    down_idx  = cell_data < 0;
        
    % Get up- and down-regulated genes
    upgenes    = ensembl(up_idx);
    downgenes  = ensembl(down_idx);
    
    % Force rho to be 10
    rho = repelem(10, length(upgenes));
    
    % Fit models and get flux distribution
    [soln, ~, ~, cell_mdl] = constrain_flux_regulation(model, ...
                                                       upgenes, downgenes, ...
                                                       kap, rho, ...
                                                       eps, ...
                                                       isgenes, ...
                                                       pfba);
    
    [geneKO, rxnKO] = knockOut(cell_mdl, 'All');
    
    % Save data
    scRxnKO(:, j)  = rxnKO;
    scGeneKO(:, j) = geneKO
    scRxnFx(:, j)  = soln;
    
    save(filename, ...
         "scRxnKO", "scGeneKO", "scRxnFx", "j", ...
         "-append"); 
end