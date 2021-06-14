%% Single-cell Flux Balance Analysis Parameter Sensitivity Analysis (PSA)
% *Author:* Scott Campit
%% Summary
% This notebook computes 5 different metabolic flux and growth rate profiles 
% by systematically modifying rho, kappa, and kappa2.
%% 
% # Default:              Kappa = 1;      Rho = 1;       
% # Lowest Kappa:  Kappa = 1E-3; Rho = 1; 
% # Lower Kappa:    Kappa = 1E-2; Rho = 1;
% # Low Kappa:       Kappa = 1E-1; Rho = 1;     
% # Lowest Rho:      Kappa = 1;       Rho = 1E-3;
% # Lower Rho:        Kappa = 1;       Rho = 1E-2;
% # Low Rho:           Kappa = 1;       Rho = 1E-1;        
%% 
% This assumes kappa2 is set to the default parameter of 1E-6. Additionally, 
% the default values for rho and kappa are 1.
%% RECON1 PSA
% This section performs parameter sensitivity analysis for RECON1. 
% A. Load the COBRA Toolbox

% Initialize metabolic modeling components
clear all;

% Needed to load COBRA onto Great lakes
addpath('/nfs/turbo/umms-csriram/scampit/Software/cobratoolbox');
initCobraToolbox(false); changeCobraSolver('gurobi', 'all');
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
% C. Set flux balance analysis hyperparameters
% Assign default values of parameters for the iMAT algorithm.

eps = [];
kap = [];
rho = [];
isgenes = true;
eps2 = [];
pfba = true;
% D. Set up parallel toolbox
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
parpool;
% E. Perform PSA
% Default

rho = 1;
kap = 1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/default/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');
% Kappa Perturbations

rho = 1;

% Kappa = 1E-3
kap = 1E-3;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/kappa1E-3/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');
% Kappa = 1E-2
kap = 1E-2;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/kappa1E-2/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');

% Kappa = 1E-1
kap = 1E-1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/kappa1E-1/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');
% Rho Perturbations

kap = 1;

% Rho = 1E-3
rho = 1E-3;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/rho1E-3/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');

% Rho = 1E-2
rho = 1E-2;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/rho1E-2/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');

% Rho = 1E-1
rho = 1E-1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/ihuman/rho1E-1/";
[~, ~, ~, ~, ~] = scfba(model, ensembl_zdata, time_course, ensembl, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'iHUMAN');
%% Summary
% This notebook goes through bulk and single cell flux balance analysis. To 
% get differentially active or differentially sensitive metabolic reactions, you 
% can check out the |scDRSA.mlx| file.