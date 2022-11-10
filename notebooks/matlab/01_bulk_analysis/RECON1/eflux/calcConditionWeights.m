function [rxns, weights] = calcConditionWeights(model, expression, params)
%calcConditionWeights calculates the weights using calculateWeight for all
%reactions within a given condition. COBRA Toolbox is a dependency.
% model:COBRA structure object. 
% expression:structure. Must contain `ids`, `pvalues` and `effect` fields. 
% params:structure. Optional. Contains pre-defined hyperparameters. 
    [~, ~, rxns] = deleteModelGenes(model, cellstr(expression.ids));
    p = expression.pvalues;
    z = expression.effect;
    
    try
        epsilon = params.epsilon;
    catch ME
        epsilon = repmat(1e-6, size(p));
    end

    weights = zeros(size(z));
    for i = 1:length(z)
        zi = z(i);
        pi = p(i);
        ep = epsilon(i);
        weights(i) = calculateWeight(zi, pi, ep);
    end
    
end