% `constrain_flux_regulation` uses the iMAT algorithm to constrain the
% metabolic model using up-regulated or down-regulated genes / reactions.
% @author: Sriram Chandrasekaran
function [constrained_model, solution] =  constrain_flux_regulation...
        (unconstrained_model, onreactions, offreactions, kappa, rho, epsilon, ...
        mode, epsilon2, minfluxflag)

if (~exist('mode','var')) || (isempty(mode))
        mode = 1;
end
if mode == 0 % genes
        [~,~,onreactions,~] =  deleteModelGenes(unconstrained_model, onreactions);
        [~,~,offreactions,~] =  deleteModelGenes(unconstrained_model, offreactions);
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
        minfluxflag = true;
end

% Parsimonious flux balance analysis
if minfluxflag
        kappa = [kappa(:); ones(size(setdiff(unconstrained_model.rxns, offreactions)))*1E-6];
        epsilon2 = [epsilon2; zeros(size(setdiff(unconstrained_model.rxns, offreactions)))];
        offreactions = [offreactions(:); setdiff(unconstrained_model.rxns, offreactions)];
end
%BIOMASS_OBJ_POS = find(ismember(unconstrained_model.rxns, 'BIOMASS_reaction'));

model = unconstrained_model;
model.A = unconstrained_model.S;
model.obj = unconstrained_model.c;
model.rhs = unconstrained_model.b;
if exist('model1.csense','var') && ~isempty(unconstrained_model.csense)
        model.sense = unconstrained_model.csense;
        model.sense(ismember(model.sense,'E')) = '=';
        model.sense(ismember(model.sense,'L')) = '<';
        model.sense(ismember(model.sense,'G')) = '>';
else
        model.sense =repmat( '=',[size(unconstrained_model.S,1),1]);
end
model.lb = unconstrained_model.lb;
model.ub = unconstrained_model.ub;
model.vtype = repmat('C', size(unconstrained_model.S, 2), 1);
model.modelsense = 'max';
nrows = size(model.A, 1);
ncols = size(model.A, 2);
M = 10000;
objpos = find(unconstrained_model.c);
number_of_rxns = length(unconstrained_model.rxns);

for on_rxns = 1:length(onreactions)
        rxnpos = find(ismember(unconstrained_model.rxns,onreactions(on_rxns)));
        
        % xi - (eps + M)ti >= -M
        %               ti = 0 or 1.
        
        rowpos = size(model.A,1) + 1;
        colpos = size(model.A,2) + 1;
        
        model.A(rowpos, rxnpos) = 1;
        model.A(rowpos, colpos) = -(1*epsilon(on_rxns) + M);
        model.rhs(rowpos) = -M;
        model.sense(rowpos) = '>';
        model.vtype(colpos) = 'B';
        model.obj(colpos) = 1*rho(on_rxns);
        model.lb(colpos) = 0;
        model.ub(colpos) = 1;
        
        % xi + (eps + M)ri <= M
        %               ri = 0 or 1.
        
        rowpos = size(model.A, 1) + 1;
        colpos = size(model.A, 2) + 1;
        
        model.A(rowpos, rxnpos) = 1;
        model.A(rowpos, colpos) = (1*epsilon(on_rxns) + M);
        model.rhs(rowpos) = M;
        model.sense(rowpos) = '<';
        model.vtype(colpos) = 'B';
        model.obj(colpos) = 1*rho(on_rxns);
        model.lb(colpos) = 0;
        model.ub(colpos) = 1;
end

% Constraints for off reactions: flux is minimized.
% Soft constraints can be violated if neccesary, and these reactions
% can carry some flux. Higher magnitude (kappa), higher penalty.

for off_rxn = 1:length(offreactions)
        rxnpos = find(ismember(unconstrained_model.rxns, offreactions(off_rxn)));
        % xi + si >= -eps2
        %      si >= 0
        % rho(ri + si)
        
        rowpos = size(model.A, 1) + 1;
        colpos = size(model.A, 2) + 1;
        model.A(rowpos, rxnpos) = 1;
        model.A(rowpos, colpos) = 1;
        model.rhs(rowpos) = -epsilon2(off_rxn);
        model.sense(rowpos) = '>';
        model.vtype(colpos) = 'C';
        model.lb(colpos) = 0;
        model.ub(colpos) = 1000;
        model.obj(colpos) = -1*kappa(off_rxn);
        
        % constraint 2
        % xi - ri <= eps2
        %      ri >= 0
        rowpos = size(model.A, 1) + 1;
        colpos = size(model.A, 2) + 1;
        model.A(rowpos, rxnpos) = 1;
        model.A(rowpos, colpos) = -1;
        model.rhs(rowpos) = epsilon2(off_rxn);
        model.sense(rowpos) = '<';
        model.vtype(colpos) = 'C';
        model.lb(colpos) = 0;
        model.ub(colpos) = 1000;
        model.obj(colpos) = -1*kappa(off_rxn);
end

params.outputflag = 0;
constrained_model = model;
solution = gurobi(constrained_model, params);
solution.flux = solution.x(1:number_of_rxns);
solution.grate = solution.x(BIOMASS_OBJ_POS);
solution.solverObj = solution.objval;

end

