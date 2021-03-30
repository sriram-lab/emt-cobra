%% RECON1 A549 COBRA Analysis
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
% The code below is for running this script on Great Lakes

% Initialize the Umich Cluster profiles
%setupUmichClusters

% We get from the environment the number of processors
%NP = str2num(getenv('SLURM_NTASKS'));

% Create the pool for parfor to use
%thePool = parpool('current', NP);
%poolobj = parpool;
%addAttachedFiles(poolobj, {})
parpool
%% 
% This codeblock is for running this script locally.

%workers = 4;
%parpool("local", workers);
% B. Initiate data structures
% We'll use the RECON1 metabolic reconstruction to perform flux balance analysis 
% and knockouts.

% ACLX
%basepath = "~/Analysis/EMT/recon1/";

% NFS
basepath = "~/Turbo/scampit/Analysis/EMT/recon1/";

% DELL
%basepath = "D:/Analysis/EMT/recon1/";
% C. Perform Single-cell simulations
% Now let's perform some knockouts. I will perform gene knockouts to save time.

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
    [soln, ~, ~, cell_mdl] = constrain_flux_regulation(model, ...
                                                       upgenes, downgenes, ...
                                                       kap, rho, ...
                                                       eps, ...
                                                       isgenes, ...
                                                       pfba);
    [geneKO, rxnKO]  = knockOut(cell_mdl, 'All');
    
    % Save data as individual files
    cobra_cell        = struct();
    cobra_cell.id     = time_course{j}
    cobra_cell.geneko = geneKO;
    cobra_cell.rxnko  = rxnKO;
    cobra_cell.flux   = soln;
    parsave(sprintf(strcat(basepath, '%d.mat'), j), 'cobra_cell', 'j');   
end
% Bulk simulation
% Now let's perform some knockouts. I will perform gene knockouts to save time.

% rxnko  = cell(length(scDE), 1);
% flux   = cell(length(scDE), 1);
% geneko = cell(length(scDE), 1);
%% 
% This is the save file path.

% % Linux
% %filename = "~/Analysis/EMT/bulk_RNASeq_recon1_profiles.mat";
% 
% % NFS
% filename = "~/Turbo/scampit/Analysis/EMT/bulk_RNASeq_recon1_profiles.mat";
% 
% % DELL
% %filename = "D:/Analysis/EMT/bulk_RNASeq_recon1_profiles.mat";
% 
% % Save it
% save(filename, 'rxnko', 'flux', 'geneko');
%% 
% This simultaneously computes metabolic fluxes, gene KO and reaction KO data.

% parfor j = 1:length(scDE)
%     
%     expt = scDE{j};
%             
%     % Get up- and down-regulated genes
%     upgenes    = expt.up;
%     downgenes  = expt.down;
%     
%     upgenes(ismissing(upgenes)) = [];
%     downgenes(ismissing(upgenes)) = [];
%     model.genes = cellstr(model.genes);
%     
%     upgenes = cellstr(upgenes);
%     downgenes = cellstr(downgenes);
%     
%     % Force rho to be 10
%     rho = repelem(10, length(upgenes));
%     
%     % Fit models and get flux distribution
%     [soln, ~, ~, cell_mdl] = constrain_flux_regulation(model, ...
%         upgenes, downgenes, ...
%         kap, rho, ...
%         eps, ...
%         isgenes, ...
%         pfba);
%     
%     [geneKO, rxnKO] = knockOut(cell_mdl, 'All');
%     
%     % Save data for each time point as a single matrix
%     flux(j, :) = soln;
%     geneko(j, :) = geneKO;
%     rxnko(j, :) = rxnKO;
%     
%     save(filename, ...
%             "flux", "geneko", "rxnko", "j", ...
%             "-append"); 
%     
% end
%% Summary
% This notebook goes through bulk and single cell flux balance analysis. To 
% get differentially active or differentially sensitive metabolic reactions, you 
% can check out the |scDRSA.mlx| file.