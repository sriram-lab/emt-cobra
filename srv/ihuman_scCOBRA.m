%% Single Cell COBRA Analysis
% *Author*: Scott Campit
%% Summary
% This script computes the reaction knockouts and metabolic fluxes for the scCOBRA 
% models.
%% 1. Load Data
% A. Load the COBRA Toolbox

% Initialize metabolic modeling components
clear all;

% Needed to load COBRA onto Great lakes
addpath('/nfs/turbo/umms-csriram/scampit/Software/cobratoolbox');

initCobraToolbox(false); 
changeCobraSolver('gurobi', 'all');
% B. Load iHUMAN reconstruction and gene expression dataset
% First, we'll load the A549 MAGIC-imputed gene expression dataset and the iHUMAN 
% metabolic reconstruction. 
% 
% For more information, see the |save_diffexp_data.mlx| livescript to see how 
% the initial dataset was preprocessed.

% DELL
%load D:/Analysis/EMT/scRNASeq_ihuman_diffexp.mat

% ACLX
%load ~/Analysis/EMT/scRNASeq_ihuman_diffexp.mat

% NFS
load ~/Turbo/scampit/Analysis/EMT/scRNASeq_ihuman_diffexp.mat
% Set the objective function to maximize biomass

% Set up objective function
model.c = zeros(size(model.c));
model.c(string(model.rxns) == 'biomass_human') = 1;
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
% A. Setting up parallel computations
% We can set the number of CPUs to use for our calculations.
% 
% The code block below is for running this script on Great Lakes

%%%%  Initialize the Umich Cluster profiles
%setupUmichClusters

%%%%  We get from the environment the number of processors
%NP = str2num(getenv('SLURM_NTASKS'));

%%%%  Create the pool for parfor to use
%thePool = parpool('current', NP);
%poolobj = parpool;
%addAttachedFiles(poolobj, {})
parpool
%% 
% The code block below is for running this script locally.

%workers = 4;
%parpool("local", workers);
% B. Initiate data structures
% We'll use the iHUMAN metabolic reconstruction to perform flux balance analysis 
% and knockouts.

% ACLX
%basepath = "~/Analysis/EMT/human1/";

% NFS
basepath = "~/Turbo/scampit/Analysis/EMT/human1/";

% DELL
%basepath = "D:/Analysis/EMT/human1/";
% C. Perform Single-cell simulations
% Now let's perform some knockouts. I will perform gene knockouts to save time.

parfor j = 1:size(ensembl_zdata, 2)
    cell_data = ensembl_zdata(:, j);
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
    
    [geneKO, ~] = knockOut(cell_mdl, 'GeneKO');

    % Save data as individual files
    cobra_cell = struct();
    cobra_cell.id     = time_course{j}
    cobra_cell.geneko = geneKO;
    cobra_cell.flux   = soln;
    save(sprintf(strcat(basepath, '%d.mat'), j), 'cobra_cell', 'j'); 
end
%% Summary
% This notebook goes through bulk and single cell flux balance analysis. To 
% get differentially active or differentially sensitive metabolic reactions, you 
% can check out the |scDRSA.mlx| file.