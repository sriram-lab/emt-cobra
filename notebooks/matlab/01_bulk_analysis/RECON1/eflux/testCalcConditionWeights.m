function testCalcConditionWeights()
%testCalcConditionWeights tests the calcConditionWeights function.

    fprintf('Test 1. Five reaction system sampled from a random distribution \n');
    
    ids = {'rxn1', 'rxn2', 'rxn3', 'rxn4', 'rxn5'};
    zscore = [1.0347, 0.7269, -2.3034, -1.2939, 2.7873];
    pvalue = [0.300809, 0.467287, 0.021279, 0.196011, 0.005315];
    epsilon = repmat(1e-3, size(zscore));

    fprintf('Reaction IDs: \n');
    disp(ids);
    fprintf('Z-score vector for condition j \n');
    disp(zscore);
    fprintf('p-value vector for condition j \n');
    disp(pvalue);
    fprintf('epsilon vector \n');
    disp(epsilon);
    
    weight = zeros(size(zscore));
    for i = 1:length(ids)
        weight(i) = calculateWeight(zscore(i), pvalue(i), epsilon(i));
    end

    disp('weights: ');
    disp(weight);

end