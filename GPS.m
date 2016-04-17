%==========================================================================
%                        The GPS Algorithm
%==========================================================================
% DATE
%        April 2015
%
% AUTHOR
%        Bryan Chong
%
% REFERENCE
%        L. Sun, L. Hong, Z. Hu. INFORMS. 2014
%==========================================================================
%
% INPUT
%        x0
%              Matrix (size = 'dim' X 'NumStartingSol') of 'NumStartingSol'
%              initial solutions to the solver
%              For R-SPLINE, NumStartingSol=1
%              Each initial solution is of size dim X 1
%        problem
%              Problem function name
%        problemseed
%              Substream index (integer >=1)
%        solverseed
%              Input seed for GPS (integer between 1 and 2147483646)
%        budget
%              Vector of size NumSoln, where NumSoln is
%              the number of solutions returned by the solver
%              for example, if budget = [500 1000] then NumSoln
%              is 2 and the solver returns best available solutions after
%              every 500 calls to the oracle
%        logfilename
%
%
% OUTPUT
%        Ancalls
%              An array (size = 'NumSoln' X 1) of buget expended
%        A
%              An array (size = 'NumSoln' X 'dim') of solutions
%              returned by solver
%        Afn
%              An array (size = 'NumSoln' X 1) of estimates of expected
%              objective function value
%        AFnVar
%              An array of variances corresponding to
%              the objective function at A
%              Equals NaN if solution is infeasible
%        AFnGrad
%              An array of gradient estimates at A; not reported
%        AFnGardCov
%              An array of gradient covariance matrices at A; not reported
%        Aconstraint
%              A vector of constraint function estimators; not applicable
%        AConstraintCov
%              An array of covariance matrices corresponding to the
%              constraint function at A; not applicable
%        AConstraintGrad
%              An array of constraint gradient estimators at A; not
%              applicable
%        AConstraintGradCov
%              An array of covariance matrices of constraint gradient
%              estimators at A; not applicable
%
%==========================================================================

%% GPS
function [Ancalls, A, Afn, AFnVar, AFnGrad, AFnGradCov, ...
    Aconstraint, AConstraintCov, AConstraintGrad, ...
    AConstraintGradCov] = GPS(x0, problem, problemseed, ...
    solverseed, budget, logfilename, minmax)
%% Unreported
AFnGrad = NaN;
AFnGradCov = NaN;
Aconstraint = NaN;
AConstraintCov = NaN;
AConstraintGrad = NaN;
AConstraintGradCov = NaN;

%% Parameters
%minM = -10^10; % Minimum function mean. Doesn't work for min problems
%minSig = 10^-5; % Minimum function variance. Doesn't work for min problems
r = 30; % Number of replications for each solution
s = 10; % Number of solutions to use at each iteration
a = 3; % Used in gamma function
b = 1.5; % Used in lambda function
%T = 20; % Number of steps used in MCC sampling
sig = 2.5; % Gaussian covariance coefficient
fastFlag = 0; % 1 if use sigmaTilde as in Proposition 2 (as suggested by author - faster but not as good), 0 to use full sigmaTilde
% Note: Unsure if bug, but not much difference when using Prop 2 sigmaTilde. Leave at value 0 for most accurate results.

%% Determined from Parameters

if min(budget) < 2*s*r
    fprintf('A budget is too small for a good quality run of GPS.');
    return;
end
iters = floor(budget/(s*r) - 1); % number of iterations of the algorithm for each solution
iterCount = 1; % which budget we are currently looking at
k = max(iters); % total number of iterations to run
fprintf('Total Iterations for GPS: %d\n',k);
NumSoln = length(budget); % number of solutions returned by solver
dim = size(x0, 2); % solution dimension

%% Initialize output
Ancalls = zeros(NumSoln,1);
A = zeros(NumSoln,dim);
Afn = zeros(NumSoln,1);
AFnVar = zeros(NumSoln,1);

%% Generate Random Variables
genSeedStream = RandStream.create('mrg32k3a', 'NumStreams', 1);
% Set the substream to the "seed"
genSeedStream.Substream = solverseed;

% Generate uniform random variables for acceptance-rejection method
OldStream = RandStream.setGlobalStream(genSeedStream); % Temporarily store old stream
genSeeds = unidrnd(100000,s*k,1);

% Restore old random number stream
RandStream.setGlobalStream(OldStream);

%% Step 0
evals = 0;
uniqSS = unique(x0, 'rows');
numUniq = size(uniqSS, 1);
n0 = r*ones(numUniq, 1); % Number of replications for each solution
G0 = zeros(numUniq, r); % Values for each starting solution X replication simulated

probHandle = str2func(problem);
for i=1:numUniq
    for j=1:r
        [fn, ~, ~, ~, ~, ~, ~, ~] = probHandle(uniqSS(i,:), 1, problemseed); problemseed = problemseed + 1; evals = evals + 1;
        G0(i,j) = minmax*fn;
    end
end
Gbar0 = mean(G0,2); % function means
var0 = var(G0,0,2); % function variances
GbarPrev = Gbar0;
varPrev = var0;
nPrev = n0;
SPrev = uniqSS;

bestVals = max(GbarPrev)*ones(1,numUniq*r); % Best values at each function evaluation

for i=1:k
    %% Step 1
    %% Step 2
    S = zeros(s,dim);
    G = zeros(s,r);
    for j = 1:s
        S(j,:) = ARS(problem, SPrev, b, GbarPrev, varPrev, nPrev, sig, a, genSeeds(s*(i-1)+j), fastFlag);
        %S(j,:) = MCCS(SPrev(mod(j,length(SPrev)+1),:), SPrev, b, GbarPrev, varPrev, nPrev, sig, VarBds, T, a, genSeeds(s*(i-1)+j));
        for m = 1:r
            [fn, ~, ~, ~, ~, ~, ~, ~] = probHandle(S(j,:), 1, problemseed); problemseed = problemseed + 1; evals = evals + 1;
            G(j,m) = minmax*fn;
        end
    end
    %% Step 3
    % Collate rows and handle GBar and nReps
    for j = 1:s
        row = S(j,:);
        repeatCheck = ismember(SPrev,row,'rows');
        if sum(repeatCheck) > 0 % Solution has been used before
            ind = find(repeatCheck);
            newGBar = (GbarPrev(ind)*nPrev(ind) + sum(G(j,:))) / (nPrev(ind) + r);
            newN = nPrev(ind) + r;
            % Update variance
            % Note: (n-1)s^2_x + xbar/n_x = sum (x^2)
            sumXsquared = (nPrev(ind)-1)*varPrev(ind) + GbarPrev(ind)/nPrev(ind);
            varPrev(ind) = (sumXsquared + sumsqr(G(j,:)) - newGBar^2)/newN;
            GbarPrev(ind) = newGBar;
            nPrev(ind) = newN;
        else % New solution
            GbarPrev = [GbarPrev; sum(G(j,:))/r];
            varPrev = [varPrev; var(G(j,:))];
            nPrev = [nPrev; r];
            SPrev = [SPrev; S(j,:)];
        end
    end
    fprintf('Completed iteration %d\n', i);
    [bestVal, indBest] = max(GbarPrev);
    bestVals = [bestVals ones(1,s*r)*bestVal];
    % Update outputs, and handle the case where two budgets result in the
    % same number of iterations
    if i == iters(iterCount) && iterCount <= NumSoln
        origIter = iterCount;
        uniqIter = 0;
        while uniqIter == 0
            if iters(iterCount) == iters(iterCount) + 1
                iterCount = iterCount + 1;
            else
                uniqIter = 1;
            end
        end
        numSameBudget = length(origIter:iterCount);
        Ancalls(origIter:iterCount) = ones(numSameBudget,1)*evals;
        bSol = SPrev(indBest,:);
        A(origIter:iterCount,:) = ones(numSameBudget,1)*bSol;
        [bVal, ~, ~, ~, ~, ~, ~, ~] = probHandle(bSol, budget(ceil(length(budget)/2)), 1);
        Afn(origIter:iterCount) = ones(numSameBudget,1)*bVal;
        AFnVar(origIter:iterCount) = ones(numSameBudget,1)*varPrev(indBest);
        iterCount = iterCount + 1;
    end
end
[bestVal, indBest] = max(GbarPrev);
bestSol = SPrev(indBest,:);
bestVal = minmax*bestVal;
plot(1:evals, minmax*bestVals);
end

% Lambda function
function [lambdaOutput] = lambda(x, S, b)
% x is the new solution value
% S is the set of all simulated solutions
% b is the lambda parameter
d = size(S,1);
lambdaOutput = zeros(d,1);

denom = 0;
for i=1:d
    numer = norm(x-S(i,:))^-b;
    denom = denom + numer;
    if isequal(x,S(i,:))
        lambdaOutput(i) = 1;
    elseif sum(ismember(S, x, 'rows')) ~= 0
        lambdaOutput(i) = 0;
    else
        lambdaOutput(i) = numer;
    end
end
lambdaOutput = lambdaOutput/denom;

% Check that lambda sums to 1
if sum(lambdaOutput) - 1 > 0.001
    fprintf('Lambda components do not sum to 1');
    return;
end
end

% Two-input gamma function
function [gamma2Output] = gamma2(x1, x2, a)
gamma2Output = exp(-a*norm(x1-x2));
end

% One-input gamma function
function [gamma1Output] = gamma1(x, S, a)
d = size(S,1);
gamma1Output = zeros(d,1);
for i=1:d
    gamma1Output(i) = gamma2(x, S(i,:),a);
end
end

% Gamma matrix
function [GammaOutput] = Gamma(S,a)
d = size(S,1);
GammaOutput = zeros(d);
for i=1:d
    for j=i:d
        GammaOutput(i,j) = gamma2(S(i,:), S(j,:), a);
        GammaOutput(j,i) = GammaOutput(i,j);
    end
end
end

% Returns conditional mean and variance
function [CondMean, CondVar] = CondMeanAndVar(x, S, b, Gbar, vars, nReps, sig, a)
d = size(S,1);
if d ~= size(Gbar,1) || d ~= size(vars,1) || d ~= length(nReps)
    fprintf('Dimension mismatch in calling GPSMeanAndVar function');
    return;
end

lamb = lambda(x,S,b);
covMat = diag(vars./nReps);

CondMean = lamb'*Gbar;
CondVar = sig^2*(1-2*lamb'*gamma1(x,S,a) + lamb'*Gamma(S,a)*lamb) + ...
    lamb'*covMat*lamb;
end

% Returns a quicker conditional mean and variance, using sigmaTilde from
% Proposition 2
function [CondMean, CondVar] = FastCondMeanAndVar(x, S, b, Gbar, vars, nReps, sig, a)
d = size(S,1);
if d ~= size(Gbar,1) || d ~= size(vars,1) || d ~= length(nReps)
    fprintf('Dimension mismatch in calling GPSMeanAndVar function');
    return;
end

lamb = lambda(x,S,b);
covMat = diag(vars./nReps);

CondMean = lamb'*Gbar;

dx = min(sqrt(sum([repmat(x,[d,1])-S].^2,2)));
hx = exp(-a*dx);
CondVar = sig^2*(1-hx)^2 + lamb'*covMat*lamb;
end

% Acceptance-Rejection Sampling
function [ARSOutput] = ARS(problem, S, b, Gbar, vars, n, sig, a, ARSseed, fastFlag)

%% Generate Random Variables
[strucSeedStream, unidStream] = RandStream.create('mrg32k3a', 'NumStreams', 2);
% Set the substream to the "seed"
strucSeedStream.Substream = ARSseed;
unidStream.Substream = ARSseed;
OldStream = RandStream.setGlobalStream(strucSeedStream); % Temporarily store old stream
c = max(Gbar);
probStrucHandle = str2func([problem 'Structure']);
strucSeed = unidrnd(10000000);
[~, ~, ~, ~, ~, ~, ~, sols, ~, ~, ~] = probStrucHandle(1000, strucSeed);
count = 1;
found = 0;

RandStream.setGlobalStream(unidStream);
while found == 0
    U = rand();
    sol = sols(count,:); count = count + 1;
    if fastFlag == 1
        [condMean, condVar] = FastCondMeanAndVar(sol, S, b, Gbar, vars, n, sig, a);
    else
        [condMean, condVar] = CondMeanAndVar(sol, S, b, Gbar, vars, n, sig, a);
    end
    if count > 1000
        strucSeed = strucSeed + 1;
        [~, ~, ~, ~, ~, ~, ~, sols, ~, ~, ~] = probStrucHandle(1000, strucSeed);
        count = 1;
    end
    Prob = 1-normcdf(c,condMean,condVar);
    if U <= 2*Prob
        found = 1;
        ARSOutput = sol;
    end
end
RandStream.setGlobalStream(OldStream);
end

%%  Alternative Sampling Algorithms
% Markov Chain Coordinate Sampling
function [MCCSOutput] = MCCS(x0, S, b, Gbar, vars, n, sig, VarBds, T, a, MCCSseed)
if sum(VarBds(:,2) - VarBds(:,1)) == inf
    fprintf('A variable is unbounded.')
    return;
end
%% Generate Random Variables
[dimStream, unifStream, ranStream] = RandStream.create('mrg32k3a', 'NumStreams', 3);
d = length(x0);
% Set the substream to the "seed"
dimStream.Substream = MCCSseed;
unifStream.Substream = MCCSseed;
ranStream.Substream = MCCSseed;

% Generate uniform random variables for acceptance-rejection method
OldStream = RandStream.setGlobalStream(dimStream); % Temporarily store old stream
Dims = unidrnd(d,T,1);

% Generate exponential random variables for acceptance-rejection method
RandStream.setGlobalStream(unifStream);
Unifs = rand(T,1);

% Generate probability of a "trip" during a failure
RandStream.setGlobalStream(ranStream);

c = max(Gbar);
y = x0;
for t = 1:T
    c1 = ceil(VarBds(Dims(t),1));
    c2 = floor(VarBds(Dims(t),2));
    ran = unidrnd(c2-c1) + c1 - 1;
    if ran > y(Dims(t))-1
        ran = min(c2, ran + 1);
    end
    z = y;
    z(Dims(t)) = ran;
    
    [yMean, yVar] = CondMeanAndVar(y, S, b, Gbar, vars, n, sig, a);
    [zMean, zVar] = CondMeanAndVar(z, S, b, Gbar, vars, n, sig, a);
    yProb = 1-normcdf(c,yMean,yVar);
    zProb = 1-normcdf(c,zMean,zVar);
    
    if Unifs(t) <= zProb/yProb
        y = z;
    end
end

% Restore old random number stream
RandStream.setGlobalStream(OldStream);

MCCSOutput = y;
end

% Generate starting solutions for AR sampling
function [StartingSols] = ARSols(NumSols, VarBds, IntVars, SolSeed)
dim = size(VarBds,1);
StartingSols = zeros(NumSols, dim);

[solStream] = RandStream.create('mrg32k3a', 'NumStreams', 1);
solStream.Substream = SolSeed;
OldStream = RandStream.setGlobalStream(solStream); % Temporarily store old stream
for i=1:dim
    if IntVars(i) == 0 % Real variables
        StartingSols(:,i) = unifrnd(VarBds(i,1), VarBds(i,2), NumSols, 1);
    elseif IntVars(i) == 1 % Integer variables
        StartingSols(:,i) = unidrnd(VarBds(i,2) - VarBds(i,1) + 1, NumSols, 1) + VarBds(i,1) - 1;
    else
        fprintf('GPS cannot handle problems with categorical variables');
        return;
    end
end
RandStream.setGlobalStream(OldStream);
end

% Acceptance-Rejection Sampling without calling structure
function [ARSOutput] = ARS2(S, b, Gbar, vars, n, sig, a, ARSseed, VarBds, IntVars)

%% Generate Random Variables
[strucSeedStream, unidStream] = RandStream.create('mrg32k3a', 'NumStreams', 2);
% Set the substream to the "seed"
strucSeedStream.Substream = ARSseed;
unidStream.Substream = ARSseed;
OldStream = RandStream.setGlobalStream(strucSeedStream); % Temporarily store old stream
c = max(Gbar);
Solseed = ARSseed;
sols = ARSols(1000, VarBds, IntVars, Solseed);
count = 1;
found = 0;

RandStream.setGlobalStream(unidStream);
while found == 0
    U = rand();
    sol = sols(count,:); count = count + 1;
    [condMean, condVar] = CondMeanAndVar(sol, S, b, Gbar, vars, n, sig, a);
    if count > 1000
        Solseed = Solseed + 1;
        sols = ARSols(1000, VarBds, IntVars, Solseed);
        count = 1;
    end
    Prob = 1-normcdf(c,condMean,condVar);
    if U <= 2*Prob
        found = 1;
        ARSOutput = sol;
    end
end
RandStream.setGlobalStream(OldStream);
end
