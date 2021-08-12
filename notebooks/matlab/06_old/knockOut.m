function [geneKO, rxnKO] = knockOut(model, KOType)
%% KNOCKOUT Knock Out 
% The |knockOut| function uses a COBRA model to compute the growth rate from 
% either gene or reaction knockouts. 
% 
% *INPUT*
% 
% |model|:                A structure of a COBRA-formatted metabolic model.
% 
% |KOType|:              A string flag specifying the knockout type. 
% 
% *OUTPUT*
% 
% |geneKOGrate|:   A vector containing the growth rates corresponding to single 
% gene knockouts
% 
% |rxnKOGrate|:     A vector containing the growth rate corresponding to single 
% reaction knockouts
    
    % GUROBI parameters in the model
    params.outputflag          = 0;
    params.Threads             = 5;
    params.Seed                = 314;
    params.NumericFocus        = 3;
        
    function geneKOGrate = singleGeneKO(model, params)
        uniqueMetabolicGenes      = unique(string(model.genes));
        parfor g = 1:length(uniqueMetabolicGenes)
            tmp                  = model;
            try
                tmp              = deleteModelGenes(tmp, char(uniqueMetabolicGenes(g)));
                solution         = gurobi(tmp, params);
                geneKOGrate(g)   = solution.x(tmp.c == 1);
            catch ME % Infeasible flux solutions
                warning(strcat('Error when knocking out gene:', string(g)));
                geneKOGrate(g) = NaN;
            end
        end
    end
    function rxnKOGrate = singleRxnKO(model, params)
        rxnKOGrate               = zeros(size(model.rxns)); 
        rxnKOFlux                = zeros(length(model.rxns), length(model.rxns));
        parfor r = 1:length(model.rxns)
            tmp                  = model;
            tmp.lb(r)            = 0; 
            tmp.ub(r)            = 0;
            try
                solution         = gurobi(tmp, params);
                rxnKOGrate(r, 1) = solution.x(tmp.c == 1);
            catch ME % Infeasible flux solutions
                warning(strcat('Error when knocking out reaction:', string(r)));
                rxnKOGrate(r, 1) = NaN;
            end
        end
    end
    switch KOType
        case 'GeneKO'
            geneKO = singleGeneKO(model, params);
            rxnKO  = [];
        case 'RxnKO'
            geneKO = [];
            rxnKO  = singleRxnKO(model, params);
        case 'All'
            geneKO = singleGeneKO(model, params);
            rxnKO  = singleRxnKO(model, params);
    end
end