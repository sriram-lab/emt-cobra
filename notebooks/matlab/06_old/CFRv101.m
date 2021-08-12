function [CFR_mdl, soln] = CFRv101(model, expression, params)
%% CFRV101 Constrained Flux Regulation (CFR) V1.01
% CFR uses expression data as a constraint to compute accurate metabolic fluxes.
% 
% *Updates from previous version (S.E.C):*
%% 
% # CFR script is broken up into subfunctions for modularity
% # Input metabolomics dataset and parameters are structured as structure/fields 
% to pass in fewer variables during CFR function call.
% # Updates both COBRA and Gurobi structures to pass into other optimization 
% functions and the COBRA toolbox.
    
    % Get on/off reactions from genes
    function [onreactions, offreactions] = mapGenes2Rxns(model, expression)
        upgenes = expression.upgenes;
        downgenes = expression.downgenes;
        [~, ~, onreactions] =  deleteModelGenes(model, upgenes);
        [~, ~, offreactions] =  deleteModelGenes(model, downgenes);
    end
    
    % Set constraints to be appropriate size
    function [epsilon, epsilon2, rho, kappa] = unpackModelConstraints(model, params, onreactions, offreactions)
        epsilon = params.epsilon;
        if (~exist('epsilon','var')) || (isempty(epsilon))
            warning("No epsilon was provided. Setting default value to 1E-3.")
            epsilon = ones(size(onreactions)) * 1E-3;
        elseif numel(epsilon) == 1 
            epsilon = repmat(epsilon, size(onreactions));
        end
        
        epsilon2 = params.epsilon2;
        if (~exist('epsilon2','var')) || (isempty(epsilon2))
            warning("No epsilon2 was provided. Setting default value to 0.")
            epsilon2 = zeros(size(offreactions));
            epsilon2 = [epsilon2; zeros(size(setdiff(model.rxns, offreactions)))];
        end
        
        rho = params.rho;
        if (~exist('rho','var')) || (isempty(rho))
            warning("No rho was provided. Setting default value to 1.")
            rho = ones(size(onreactions));
        elseif numel(rho) == 1 
            rho = repmat(rho, size(onreactions));
        end
        
        kappa = params.kappa;
        if (~exist('kappa','var')) || (isempty(kappa))
            warning("No kappa was provided. Setting default value to 1.")
            kappa = ones(size(offreactions));
        elseif numel(kappa) == 1 
            kappa = repmat(kappa, size(onreactions));
        end
    end
    
    % Turn on pFBA
    function [offreactions, kappa] = setUpPFBA(model, offreactions, kappa, params)
        kappa2 = params.kappa2;
        if (~exist('kappa2','var')) || (isempty(kappa2))
            warning("No kappa2 was provided. Setting default value to 1E-6.")
            kappa = [kappa(:); ones(size(setdiff(model.rxns, offreactions))) * 1E-6];
        elseif kappa2 == 0
            warning("Kappa2 is set to 0. No pFBA.");
            kappa = [kappa(:); zeros(size(setdiff(model.rxns, offreactions)))];
        elseif numel(kappa2) == 1 
            kappa = [kappa(:); ones(size(setdiff(model.rxns, offreactions))) * kappa2];
        end
        offreactions = [offreactions(:); setdiff(model.rxns, offreactions)];
    end
    
    % Add Gurobi variables
    function gurobi_model = makeGurobiObj(model)
        if ~isfield('model', 'A')
            tmpModel            = model;
            tmpModel.A          = tmpModel.S;
            tmpModel.obj        = tmpModel.c;
            tmpModel.rhs        = tmpModel.b;
            if ~isfield('model', 'csense')
                tmpModel.sense      = repmat('=', [size(tmpModel.S, 1), 1]);
                tmpModel.csense      = tmpModel.sense;
            else
                tmpModel.sense = tmpModel.csense;
                tmpModel.sense(ismember(tmpModel.sense,'E')) = '=';
                tmpModel.sense(ismember(tmpModel.sense,'L')) = '<';
                tmpModel.sense(ismember(tmpModel.sense,'G')) = '>';
            end
            tmpModel.lb         = tmpModel.lb;
            tmpModel.ub         = tmpModel.ub;
            tmpModel.vtype      = repmat('C', [size(tmpModel.S, 2), 1]);
            tmpModel.modelsense = 'max';
            gurobi_model        = tmpModel;
        else
            gurobi_model = model;
        end
    end
    
    % Add onreaction constraints
    function tmpMdl1 = addUpregulationConstraints(gurobi_model, onreactions, rho, epsilon)
        M = 10000;
        tmpMdl1 = gurobi_model;
        
        for j = 1:length(onreactions)
    
            rxnpos = find(ismember(tmpMdl1.rxns, onreactions(j)));
            rowpos = size(tmpMdl1.A, 1) + 1;
            colpos = size(tmpMdl1.A, 2) + 1;
            
            tmpMdl1.A(rowpos, rxnpos) = 1;
            tmpMdl1.A(rowpos, colpos) = -1 * (epsilon(j) + M);
            tmpMdl1.S(rowpos, rxnpos) = 1;
            tmpMdl1.S(rowpos, colpos) = -1 * (epsilon(j) + M);
            tmpMdl1.rhs(rowpos)       = -M;
            tmpMdl1.b(rowpos)         = -M;
            tmpMdl1.sense(rowpos)     = '>';
            tmpMdl1.vtype(colpos)     = 'B';
            tmpMdl1.obj(colpos)       = rho(j);
            tmpMdl1.c(colpos)         = rho(j);
            tmpMdl1.lb(colpos)        = 0;
            tmpMdl1.ub(colpos)        = 1;
            
            rowpos = size(tmpMdl1.A, 1) + 1;
            colpos = size(tmpMdl1.A, 2) + 1;
            
            tmpMdl1.A(rowpos, rxnpos)  = 1;
            tmpMdl1.A(rowpos, colpos)  = epsilon(j) + M;
            tmpMdl1.S(rowpos, rxnpos)  = 1;
            tmpMdl1.S(rowpos, colpos)  = epsilon(j) + M;
            tmpMdl1.rhs(rowpos)        = M;
            tmpMdl1.b(rowpos)          = M;
            tmpMdl1.sense(rowpos)      = '<';
            tmpMdl1.vtype(colpos)      = 'B';
            tmpMdl1.obj(colpos)        = rho(j);
            tmpMdl1.c(colpos)          = rho(j);
            tmpMdl1.lb(colpos)         = 0;
            tmpMdl1.ub(colpos)         = 1;
            
        end
    end
    
    % Add offreaction constraints
    function tmpMdl2 = addDownregulationConstraints(tmpMdl1, offreactions, kappa, epsilon2)
        tmpMdl2 = tmpMdl1;
        for jj = 1:length(offreactions)
            rxnpos = find(ismember(tmpMdl2.rxns, offreactions(jj)));
            rowpos = size(tmpMdl2.A, 1) + 1;
            colpos = size(tmpMdl2.A, 2) + 1;
            tmpMdl2.A(rowpos, rxnpos) = 1;
            tmpMdl2.A(rowpos, colpos) = 1;
            tmpMdl2.S(rowpos, rxnpos) = 1;
            tmpMdl2.S(rowpos, colpos) = 1;
            tmpMdl2.rhs(rowpos)       = -1 * epsilon2(jj);
            tmpMdl2.b(rowpos)         = -1 * epsilon2(jj);
            tmpMdl2.sense(rowpos)     = '>';
            tmpMdl2.vtype(colpos)     = 'C';
            tmpMdl2.lb(colpos)        = 0;
            tmpMdl2.ub(colpos)        = 1000;
            tmpMdl2.obj(colpos)       = -1 * kappa(jj);
            tmpMdl2.c(colpos)         = -1 * kappa(jj); 
            
            rowpos = size(tmpMdl2.A, 1) + 1;
            colpos = size(tmpMdl2.A, 2) + 1;
            tmpMdl2.A(rowpos, rxnpos) = 1;
            tmpMdl2.A(rowpos, colpos) = -1;
            tmpMdl2.S(rowpos, rxnpos) = 1;
            tmpMdl2.S(rowpos, colpos) = -1;
            tmpMdl2.rhs(rowpos)       = epsilon2(jj);
            tmpMdl2.b(rowpos)         = epsilon2(jj);
            tmpMdl2.sense(rowpos)     = '<';
            tmpMdl2.vtype(colpos)     = 'C';
            tmpMdl2.lb(colpos)        = 0;
            tmpMdl2.ub(colpos)        = 1000;
            tmpMdl2.obj(colpos)       = -1 * kappa(jj); 
            tmpMdl2.c(colpos)         = -1 * kappa(jj); 
        end
    end
    
    % Solve using Gurobi Optimizer
    function soln = solveConstrainedFluxes(model, params2)
        soln = gurobi(model, params2);
    end
    
    % Run CFR functions
    [onreactions, offreactions]     = mapGenes2Rxns(model, expression);
    [epsilon, epsilon2, rho, kappa] = unpackModelConstraints(model, params, ...
                                                             onreactions, offreactions);
    [offreactions, kappa]           = setUpPFBA(model, offreactions, ...
                                                kappa, params);
    gurobi_model                    = makeGurobiObj(model);
    tmpMdl1                         = addUpregulationConstraints(gurobi_model, onreactions, ...
                                                                 rho, epsilon);
    CFR_mdl                         = addDownregulationConstraints(tmpMdl1, offreactions, ...
                                                                   kappa, epsilon2);
    % Ensure no summary output from Gurobi.
    params2                         = struct();
    params2.OutputFlag              = false;
    soln                            = solveConstrainedFluxes(CFR_mdl, params2);
end