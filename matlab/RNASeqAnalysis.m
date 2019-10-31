load recon1.mat

array = readcell('./../pyscripts/GSE17708mapped.xlsx');
vals = cell2mat(array(2:end, 2:end));
geneNames = array(2:end, 1);
conditions = array(1, 2:end);

muCntrl = mean(vals(:, 1:3), 2);
mu05 = mean(vals(:, 4:6), 2);
mu1 = mean(vals(:, 7:9), 2);
mu2 = mean(vals(:, 10:11), 2);
mu4 = mean(vals(:, 12:14), 2);
mu8 = mean(vals(:, 15:17), 2);
mu16 = mean(vals(:, 18:20), 2);
mu24 = mean(vals(:, 21:23), 2);
mu72 = mean(vals(:, 24:26), 2);

% difference between control and time points
muArray = [mu05, mu1, mu2, mu4, mu8, mu16, mu24, mu72];
normMuArray = muArray ./ muCntrl;
[U, ~, idx] = unique(string(geneNames), 'stable');

% Average duplicate values
for i = 1:size(normMuArray,2)
    averageMuArray(:, i) = accumarray(idx, normMuArray(:,i), [], @mean);
end

% Differentially expressed genes with respect to control
for condition = 1:size(averageMuArray, 2)
    mu = mean(averageMuArray(:, condition));
    sigma = std(averageMuArray(:, condition));
    
    for gene = 1:size(averageMuArray, 1)
        [~, pvalue(gene, condition), ~, zscore(gene, condition)]...
            = ztest(averageMuArray(gene, condition), mu, sigma);
    end
end

% Remove entries with NaN
remove = ismissing(string(geneNames));
geneNames = geneNames(~remove);
vals = vals(~remove,:);
pvalue = pvalue(~remove,:);
zscore = zscore(~remove,:);

% Get gene symbols
mask = pvalue <= 0.05;
filteredZscores = zscore;
filteredZscores(~mask) = NaN;

Tx = {'TGF05', 'TGF1', 'TGF2', 'TGF4', 'TGF8', 'TGF16', 'TGF24', 'TGF72'};
for condition = 1:size(filteredZscores, 2)    
    TxName = Tx(condition);
    
    DownPos = (filteredZscores(:, condition) <= 0);
    UpPos = (filteredZscores(:, condition) >= 0);

    UpGenes = geneNames(UpPos);
    DownGenes = geneNames(DownPos);
    
    model_ON = unique(model.genes(ismember(string(model.genes), string(UpGenes))));
    model_OFF = unique(model.genes(ismember(string(model.genes), string(DownGenes))));
    
    ON_fieldname = string(strcat(TxName, '_ON'));
    OFF_fieldname = string(strcat(TxName, '_OFF'));
    
    diffExpGenes.(ON_fieldname) = model_ON;
    diffExpGenes.(OFF_fieldname) = model_OFF;
end

% Run constrain flux regulation for the 8 models 
posToTurnOff = find(ismember(model.rxns, {'ALR', 'MGSA', 'MGSA2'}));
model.ub(posToTurnOff) = 0;
model.lb(posToTurnOff) = 0;
for condition = 1:size(Tx, 2)
    TxName = Tx(condition);
    ON_fieldname = string(strcat(TxName, '_ON'));
    OFF_fieldname = string(strcat(TxName, '_OFF'));
    
    str = strcat("[~, ", string(TxName), "Flux, ", string(TxName), "grate, ~]", ...
        "=  constrain_flux_regulation(model, diffExpGenes.(ON_fieldname), ", ...
        "diffExpGenes.(OFF_fieldname), [], [], [], 1, [], 1);");
    eval(str)
end

% Constract array
FluxArray = [TGF05Flux, TGF1Flux, TGF2Flux, TGF4Flux, ...
             TGF8Flux, TGF16Flux, TGF24Flux, TGF72Flux];
         
for condition = 1:size(FluxArray, 2)
    mu = mean(FluxArray(:, condition));
    sigma = std(FluxArray(:, condition));
    
    for reaction = 1:size(FluxArray, 1)
        [~, pvalueFlux(reaction, condition), ~, zscoreFlux(reaction, condition)]...
            = ztest(FluxArray(reaction, condition), mu, sigma);
    end
end

for condition = 1:size(FluxArray, 2)
   rxn = string(model.rxns);
   subsystem = string(model.subSystems);
   name = string(model.rxnNames);
   flux = FluxArray(:, condition);
   pval = pvalueFlux(:, condition);
   zscr = zscoreFlux(:, condition);
   array = [rxn, subsystem, name, flux, pval, zscr];
   labeledArray = ["Reaction", "Metabolic_Subsystem", "Name", "Flux", "pValue", "Zscore";
       array];
   writematrix(labeledArray, 'CFR_TS.xlsx', 'Sheet', string(Tx(condition)));
end








