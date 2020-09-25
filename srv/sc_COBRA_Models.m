%% Single Cell COBRA Analysis
% *Author*: Scott Campit
%% Summary
% This script computes the reaction knockouts and metabolic fluxes for the scCOBRA 
% models.
%% 1. Load Bulk scRNASeq Data
% This module loads the single cell differentially expressed genes |scDE.mat| 
% and the metabolic model |recon1.mat|.

% Initialize metabolic modeling components
clear all;
load ~/Analysis/EMT/scDE.mat
load ~/Data/Reconstructions/RECON1/recon1.mat

initCobraToolbox; changeCobraSolver('gurobi', 'all');
%% 2. Set up metabolic model constraints and hyperparameters
% Metabolic model constraints
% First, let's set the genes for the model to be Entrez IDs.

model = recon1;
model.genes = string(model.geneEntrezID);
%% 
% Set the objective function to maximize biomass.

% Set up objective function
model.c = zeros(size(model.c));
model.c(string(model.rxns) == 'biomass_objective') = 1;
%% 
% Let's turn off some specific reactions that cause trouble with metabolic models.

% Turn off specific reactions
posToTurnOff           = find(ismember(model.rxns, {'ALR', 'MGSA', 'MGSA2'}));
model.ub(posToTurnOff) = 0; model.lb(posToTurnOff) = 0;
% Metabolic model hyperparameters
% Assign default values of parameters for the iMAT algorithm.

% Params for constrain flux regulation. Set all to empty.
hyperparams.eps  = [];      hyperparams.isgenes = true;
hyperparams.kap  = 1;       hyperparams.eps2 = [];
hyperparams.rho  = 5;       hyperparams.pfba = true;
hyperparams.kap2 = 1E-6;    hyperparams.gamma = [];
hyperparams.hscore = false;
%% 
% We can set the number of CPUs to use for our calculations.

%workers = 4;
parpool("local", workers);
%% 
% Set file name for bulk RNASeq knockouts
% 
% scRxnKO = cell(length(scDE), 1);
% scRxnFx = cell(length(scDE), 1);
% filename = "~/Analysis/EMT/scRNASeq_bulk_profiles.mat";
% save(filename, 'scRxnKO', 'scRxnFx');
% %% 3. Compute knockout growth rates and metabolic fluxes for the aggregated scRNASeq data
% % Now let's compute the knockouts and metabolic fluxes.
% 
% for j = 1:length(scDE)
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
%     % Fit models and get 
%     [mdl, soln] = CFR(model, hyperparams, upgenes, downgenes);
%     [~, rxnKO]  = knockOut(mdl, 'RxnKO');
%     
%     % Save data
%     scRxnKO(:, j) = rxnKO;
%     scRxnFx(:, j) = soln.x(1:3744);
%     
%     save(filename, ...
%             "scRxnKO", "scRxnFx", "j", ...
%             "-append");
% end
%% 4. Load individual scRNASeq Data and preprocess

load /home/scampit/Analysis/EMT/scDE_indvidual.mat
%% 
% Read in the labels

timepoints = readcell('~/Analysis/EMT/timepoints.txt');
timepoints = string(timepoints);
%% 
% Find the cells to remove

rm_cells = find(contains(timepoints, 'rm'));
timepoints(rm_cells) = [];
%% 
% Remove the expression data in the scDE

for i = 1:length(rm_cells)
    scDE{rm_cells(i)} = [];
end
scDE = scDE(~cellfun('isempty',scDE));
%% 5. Select 100 random cells to perform knockout and fluxes

idx = datasample(1:length(scDE), 100, 'Replace', false);
for i = 1:length(idx)
    subDE{i} = scDE{idx(i)};
end
sub_timepoints = timepoints(idx);
%% 6. Compute knockout growth rates and metabolic fluxes for the aggregated scRNASeq data
% First, we'll save some empty cell arrays that will store the data.

scRxnKO = zeros([length(model.rxns), length(scDE)]);
scRxnFx = zeros([length(model.rxns), length(scDE)]);
filename = "~/Analysis/EMT/scRNASeq_indv_profiles.mat";
save(filename, 'scRxnKO', 'scRxnFx', 'rm_cells', 'sub_timepoints');
%% 
% Then, let's compute it all.

for j = 1:length(idx)
    expt = scDE{j};
            
    % Get up- and down-regulated genes
    upgenes    = expt.up;
    downgenes  = expt.down;
    
    upgenes(ismissing(upgenes)) = [];
    downgenes(ismissing(upgenes)) = [];
    model.genes = cellstr(model.genes);
    
    upgenes = cellstr(upgenes);
    downgenes = cellstr(downgenes);
    
    % Fit models and get fluxes or growth rates from knockout
    try
        [mdl, soln] = CFR(model, hyperparams, upgenes, downgenes);
        scRxnFx(:, j) = soln.x(1:3744);
    catch ME
        scRxnFx(:, j) = NaN(1:3744);
    end
    
    try
        [~, rxnKO]  = knockOut(mdl, 'RxnKO');
        scRxnKO(:, j) = rxnKO;
    catch ME
        scRxnKO(:, j) = NaN(1:3744);
    end
    
    save(filename, ...
            "scRxnKO", "scRxnFx", "j", ...
            "-append");
end