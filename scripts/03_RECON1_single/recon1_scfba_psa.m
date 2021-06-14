%% Single-cell Flux Balance Analysis Parameter Sensitivity Analysis
% *Author:* Scott Campit
%% Summary
% This notebook computes 5 different metabolic flux and growth rate profiles 
% by systematically modifying rho, kappa, and kappa2.
%% 
% # Default:         Kappa = 1;      Rho = 1;       
% # Low Kappa:  Kappa = 1E-6; Rho = 1; 
% # High Kappa: Kappa = 10;     Rho = 1     
% # Low Rho:      Kappa = 1;       Rho = 1E-6;
% # High Rho:     Kappa = 1;       Rho = 10;        
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
% B. Load RECON1 reconstruction and gene expression dataset
% This loads the data for RECON1.

% DELL
%load D:/Analysis/EMT/scRNASeq_recon1_diffexp.mat

% ACLX
%load ~/Analysis/EMT/scRNASeq_recon1_diffexp.mat

% NFS
load ~/Turbo/scampit/Analysis/EMT/scRNASeq_recon1_diffexp.mat
%% 
% Set the objective function to maximize biomass.

% Set up objective function
model.c = zeros(size(model.c));
model.c(string(model.rxns) == 'biomass_objective') = 1;
%% 
% Let's turn off some specific reactions that cause trouble with metabolic models.

%Turn off specific reactions
posToTurnOff           = find(ismember(model.rxns, {'ALR', 'MGSA', 'MGSA2'}));
model.ub(posToTurnOff) = 0; 
model.lb(posToTurnOff) = 0;
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
%parpool("local", workers);
parpool;
% E. Perform PSA
% Default

rho = 1;
kap = 1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/default/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');
% Kappa Perturbation

rho = 1;

% Kappa 1E-3
kap = 1E-3;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/kappa1E-3/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');

% Kappa 1E-2
kap = 1E-2;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/kappa1E-2/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');

% Kappa 1E-1
kap = 1E-1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/kappa1E-1/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');
% Rho Perturbation

kap = 1;

% Rho = 1E-3
rho = 1E-3;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/rho1E-3/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');
% Rho = 1E-2
rho = 1E-2;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/rho1E-2/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');
% Rho = 1E-1
rho = 1E-1;
basepath = "~/Turbo/scampit/Analysis/EMT/psa/recon1/rho1E-1/";
[~, ~, ~, ~, ~] = scfba(model, entrez_zdata, time_course, entrez, ...
                        rho, kap, eps, isgenes, pfba, ...
                        basepath, 'RECON1');
%% Summary