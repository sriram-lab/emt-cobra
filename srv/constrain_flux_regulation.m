function [fluxstate_gurobi,grate, solverobj, model]...
    =  constrain_flux_regulation(model1,onreactions,offreactions,kappa,...
    rho,epsilon,mode,minfluxflag)

if (~exist('mode','var')) || (isempty(mode))
    mode = 1;
end

if mode == 0 % genes
    [~,~,onreactions,~] =  deleteModelGenes(model1, onreactions);
    [~,~,offreactions,~] =  deleteModelGenes(model1, offreactions);
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


%if (~exist('epsilon2','var')) || (isempty(epsilon2))
%    epsilon2 = zeros(size(offreactions));
%end

if (~exist('minfluxflag','var')) || (isempty(minfluxflag))
    minfluxflag = true; % by default sum of flux through all ractions is minimized
end

params.outputflag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if minfluxflag
    kappa = [kappa(:); ones(size(setdiff(model1.rxns, offreactions)))*1E-6]; % minimize flux through all the reactions. PFBA
    %epsilon2 = [epsilon2; zeros(size(setdiff(model1.rxns, offreactions)))];
    offreactions = [offreactions(:); setdiff(model1.rxns, offreactions)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to gurobi format.
model     = model1;
model.A   = model1.S;
model.obj = model1.c;
model.rhs = model1.b;

if exist('model1.csense','var') && ~isempty(model1.csense)
    model.sense = model1.csense;
    model.sense(ismember(model.sense,'E')) = '=';
    model.sense(ismember(model.sense,'L')) = '<';
    model.sense(ismember(model.sense,'G')) = '>';
else
    model.sense =repmat( '=',[size(model1.S,1),1]);
end

model.lb = model1.lb;
model.ub = model1.ub;
model.vtype = repmat('C',size(model1.S,2),1);
model.modelsense = 'max';
nrows = size(model.A,1);
ncols = size(model.A,2);
M = 10000;
%epsilon2 = 0;
objpos = find(model1.c);
nrxns = length(model1.rxns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximize the number of reactions with proteomic or transcriptomic evidence that are ON/up-regulated.
for j = 1:length(onreactions)
    rxnpos = find(ismember(model1.rxns,onreactions(j)));
    
    %         xi - (eps + M)ti >= -M
    % ti = 0 or 1.
    
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = -(1*epsilon(j) + M);
    model.rhs(rowpos) = -M;
    model.sense(rowpos) = '>';
    model.vtype(colpos) = 'B';
    model.obj(colpos) = 1*rho(j);
    model.lb(colpos) = 0;
    model.ub(colpos) = 1;
    
    % xi + (eps + M)ri <= M
    % ri = 0 or 1.
    
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = (1*epsilon(j) + M);
    model.rhs(rowpos) = M;
    model.sense(rowpos) = '<';
    model.vtype(colpos) = 'B';
    model.obj(colpos) = 1*rho(j);
    model.lb(colpos) = 0;
    model.ub(colpos) = 1;
end

% constraints for off reactions. their flux is minimized.
% soft constraints - can be violated if neccesary - i.e some reactions
% can have some flux.  higher magnitude higher penalty

for jj = 1:length(offreactions)
    rxnpos = find(ismember(model1.rxns,offreactions(jj)));
    %        xi + si >= -eps2
    %     si >= 0
    %     rho(ri + si)
    % constraint 1
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = 1;
    model.rhs(rowpos) = 0;
    model.sense(rowpos) = '>';
    model.vtype(colpos) = 'C';
    
    % set si to be positive
    model.lb(colpos) = 0;
    model.ub(colpos) = 1000;
    model.obj(colpos) = -1*kappa(jj); % minimized
    
    % constraint 2
    %     xi - ri <= eps2
    %     ri >= 0
    % new row and column
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = -1;
    model.rhs(rowpos) = 0;
    model.sense(rowpos) = '<';
    model.vtype(colpos) = 'C';
    
    % set ri to be positive
    model.lb(colpos) = 0;
    model.ub(colpos) = 1000;
    model.obj(colpos) = -1*kappa(jj); % minimized
end

% Set params for Gurobi
params.outputflag          = 0;
params.Threads             = 2;
params.Seed                = 314;
params.NumericFocus        = 3;
solg1 = gurobi(model, params);

try
    fluxstate_gurobi = solg1.x(1:nrxns);
    grate = solg1.x(objpos);
    solverobj = solg1.objval;
    
catch
    fluxstate_gurobi = NaN;
    grate = NaN;
    solverobj = NaN;
end

end

