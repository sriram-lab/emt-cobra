function testCalculateWeight()
%testCalculateWeight tests the calculateWeight function.

    fprintf('Objectives: \n');
    fprintf('lim Z -> inf: weight -> inf \n');
    fprintf('lim Z -> 0: weight -> 1 \n');
    fprintf('lim Z -> -inf: weight -> epsilon \n\n');

    fprintf('Test 1. High positive Z-score condition \n');
    fprintf('Z-score = 10, p-value = 1e-5, epsilon = 1e-3 \n')
    zscore = 10;
    pvalue = 1e-5;
    epsilon = 1e-3;
    weight = calculateWeight(zscore, pvalue, epsilon);
    
    formatSpec = 'weight is %4.2f\n\n';
    fprintf(formatSpec, weight);

    fprintf('Test 2. High negative Z-score condition \n');
    fprintf('Z-score = -10, p-value = 1e-5, epsilon = 1e-3 \n')
    zscore = -10;
    pvalue = 1e-5;
    epsilon = 1e-3;
    weight = calculateWeight(zscore, pvalue, epsilon);
    
    formatSpec = 'weight is %4.2f\n\n';
    fprintf(formatSpec, weight);

    fprintf('Test 3. Non-significant Z-score condition \n');
    fprintf('Z-score = 1, p-value = 0.317311, epsilon = 1e-3 \n')
    zscore = 1;
    pvalue = 0.317311;
    epsilon = 1e-3;
    weight = calculateWeight(zscore, pvalue, epsilon);
    
    formatSpec = 'weight is %4.2f\n\n';
    fprintf(formatSpec, weight);

    fprintf('Test 4. Z-score close to 0 condition \n');
    fprintf('Z-score = 0, p-value = 1, epsilon = 1e-3 \n')
    zscore = 0;
    pvalue = 1;
    epsilon = 1e-3;
    weight = calculateWeight(zscore, pvalue, epsilon);
    
    formatSpec = 'weight is %4.2f\n\n';
    fprintf(formatSpec, weight);

end