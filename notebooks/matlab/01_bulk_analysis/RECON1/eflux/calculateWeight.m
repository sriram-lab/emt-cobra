function weight = calculateWeight(zscore, pvalue, epsilon)
%calculateWeight uses the Z-score and p-value to calculate a
%constraint-value for an individual reaction. epsilon is an additional
%hyperparameter that allows the constraint to remain active across all
%conditions.

% lim Z -> inf: weight -> inf
% lim Z -> 0: weight -> 1
% lim Z -> -inf: weight -> epsilon

    if zscore >= 0 
        weight = zscore * (1 - pvalue);
    else       
        weight = (1 / abs(zscore)) * (1 - pvalue);
    end

    if pvalue > 0.05
        weight = weight + 1;
    end
    
    weight = weight + epsilon; 
end