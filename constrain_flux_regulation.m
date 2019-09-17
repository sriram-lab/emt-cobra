function [constrain_model, metabolic_flux, growth_rate, solverobj]...
    =  constrain_flux_regulation(model, onreactions, offreactions, kappa,...
    rho, epsilon, mode, epsilon2, minfluxflag)

if (~exist('mode','var')) || (isempty(mode))
    mode = 1;
elseif mode == 0
    [~,~,onreactions,~] =  deleteModelGenes(model, onreactions);
    [~,~,offreactions,~] =  deleteModelGenes(model, offreactions);
end
if (~exist('epsilon','var')) || (isempty(epsilon))
    epsilon = ones(size(onreactions))*1E-3;
end
if numel(epsilon) == 1
    epsilon = repmat(epsilon, size(onreactions));
end
if (~exist('rho','var')) || (isempty(rho)) 
    rho = repmat(1, size(onreactions));
end
if numel(rho) == 1
    rho  = repmat(rho, size(onreactions));
end
if (~exist('kappa','var')) || (isempty(kappa))
    kappa = repmat(1, size(offreactions));
end
if numel(kappa) == 1
    kappa  = repmat(kappa, size(offreactions));  
end
if (~exist('epsilon2','var')) || (isempty(epsilon2))
    epsilon2 = zeros(size(offreactions));
end
if (~exist('minfluxflag','var')) || (isempty(minfluxflag))
    minfluxflag = true; % by default sum of flux through all ractions is minimized
end

% Parsimonious FBA
params.outputflag = 0;
if minfluxflag
    kappa = [kappa(:); ones(size(setdiff(model.rxns, offreactions)))*1E-6]; 
    epsilon2 = [epsilon2; zeros(size(setdiff(model.rxns, offreactions)))];
    offreactions = [offreactions(:); setdiff(model.rxns, offreactions)];
end

% Convert metabolic model to Gurobi format
gurobi_model = model;
gurobi_model.A = model.S;
gurobi_model.obj = model.c;
gurobi_model.rhs = model.b;

if exist('model1.csense','var') && ~isempty(model.csense)
    gurobi_model.sense = model.csense;
    gurobi_model.sense(ismember(model.sense,'E')) = '=';
    gurobi_model.sense(ismember(model.sense,'L')) = '<';
    gurobi_model.sense(ismember(model.sense,'G')) = '>';
else
    gurobi_model.sense = repmat( '=', [size(gurobi_model.S, 1), 1]);
end

gurobi_model.lb = model.lb;
gurobi_model.ub = model.ub;
gurobi_model.vtype = repmat('C', size(model.S, 2), 1);
gurobi_model.modelsense = 'max';
nrows = size(gurobi_model.A,1);
ncols = size(gurobi_model.A,2);
M = 10000;
objpos = find(model.c);
number_of_reactions = length(model.rxns);

% maximize the number of reactions with proteomic or transcriptomic evidence that are ON/up-regulated.
for on_rxn = 1:length(onreactions)
    rxnpos = find(ismember(gurobi_model.rxns, onreactions(on_rxn)));
    
    % xi - (eps + M)ti >= -M
    % ti = 0 or 1.
    
    new_row = size(gurobi_model.A, 1) + 1;
    new_col = size(gurobi_model.A, 2) + 1;
    
    gurobi_model.A(new_row, rxnpos) = 1;
    gurobi_model.A(new_row, new_col) = -1*(1*epsilon(on_rxn) + M);
    gurobi_model.rhs(new_row) = -1*M;
    gurobi_model.sense(new_row) = '>';
    gurobi_model.vtype(new_col) = 'B';
    gurobi_model.obj(new_col) = 1*rho(on_rxn);
    gurobi_model.lb(new_col) = 0;
    gurobi_model.ub(new_col) = 1;
    
    % xi + (eps + M)ri <= M
    % ri = 0 or 1.
    
    new_row = size(gurobi_model.A, 1) + 1;
    new_col = size(gurobi_model.A, 2) + 1;
    
    gurobi_model.A(new_row, rxnpos) = 1;
    gurobi_model.A(new_row, new_col) = (1*epsilon(on_rxn) + M);
    gurobi_model.rhs(new_row) = M;
    gurobi_model.sense(new_row) = '<';
    gurobi_model.vtype(new_col) = 'B';
    gurobi_model.obj(new_col) = 1*rho(on_rxn);
    gurobi_model.lb(new_col) = 0;
    gurobi_model.ub(new_col) = 1;
end

% constraints for off reactions. their flux is minimized.
% soft constraints - can be violated if neccesary - i.e some reactions
% can have some flux.  higher magnitude higher penalty

for off_rxn = 1:length(offreactions)
    rxnpos = find(ismember(gurobi_model.rxns, offreactions(off_rxn)));
    % xi + si >= -eps2
    % si >= 0
    % rho(ri + si)
    
    % constraint 1
    new_row = size(gurobi_model.A, 1) + 1;
    new_col = size(gurobi_model.A, 2) + 1;
    gurobi_model.A(new_row, rxnpos) = 1;
    gurobi_model.A(new_row, new_col) = 1;
    gurobi_model.rhs(new_row) = -1*epsilon2(off_rxn);
    gurobi_model.sense(new_row) = '>';
    gurobi_model.vtype(new_col) = 'C';
    
    % set si to be positive
    gurobi_model.lb(new_col) = 0;
    gurobi_model.ub(new_col) = 1000;
    gurobi_model.obj(new_col) = -1*kappa(off_rxn); % minimized
    
    % constraint 2
    % xi - ri <= eps2
    % ri >= 0

    new_row = size(gurobi_model.A,1) + 1;
    new_col = size(gurobi_model.A,2) + 1;
    gurobi_model.A(new_row, rxnpos) = 1;
    gurobi_model.A(new_row, new_col) = -1;
    gurobi_model.rhs(new_row) = epsilon2(off_rxn);
    gurobi_model.sense(new_row) = '<';
    gurobi_model.vtype(new_col) = 'C';
    
    % set ri to be positive
    gurobi_model.lb(new_col) = 0;
    gurobi_model.ub(new_col) = 1000;
    gurobi_model.obj(new_col) = -1*kappa(off_rxn); % minimized
end

constrain_model = gurobi_model;

solution = gurobi(constrain_model, params);

try
    metabolic_flux = solution.x(1:length(model.rxns));
    growth_rate = solution.x(objpos);
    solverobj = solution.objval;
catch
    metabolic_flux = NaN;
    growth_rate = NaN;
    solverobj = NaN;
end

end

