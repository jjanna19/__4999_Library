%==========================================================================
%                The Random Search Algorithm
%==========================================================================
% DATE
%        April 2016
%
% AUTHOR
%        Arielle Anderer and Anna Dong, Credited to Anna Dong's More 
%           Conservative Random Search Algorithm
%
%==========================================================================
%
% INPUT
%        x0
%              Matrix (size = 'NumStartingSol' X 'dim') of 'NumStartingSol'
%              initial solutions to the solver
%              For Random Search, NumStartingSol<=9
%              Each initial solution is of size dim X 1
%        problem
%              Problem function name
%        problemseed
%              Substream index (integer >=1)
%        solverseed
%              Input seed for Random Search (integer between 1 and 2147483646)
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
%              An array (size = 'NumSoln' X 1) of budget expended
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

%% Random Search
function [Ancalls, A, Afn, AFnVar, AFnGrad, AFnGradCov, ...
    Aconstraint, AConstraintCov, AConstraintGrad, ...
    AConstraintGradCov] = RandSearch(x0, problem, problemseed, ...
    solverseed, budget, logfilename, minmax)
%% Unreported
AFnGrad = NaN;
AFnGradCov = NaN;
Aconstraint = NaN;
AConstraintCov = NaN;
AConstraintGrad = NaN;
AConstraintGradCov = NaN;

%%
r = 30; % Runlength time for each solution
s = 10; % Number of solutions to replicate/simulate at each iteration
numGen = 1000; % Number of starting solutions to (randomly) search among at each iteration

if min(budget) < s*r
    fprintf('A budget is too small for a good quality run of RandomSearch.');
    return;
end
iters = floor(budget/(s*r)); % number of iterations of the algorithm for each solution
iterCount = 1; % which budget we are currently looking at
Ancalls = (iters*s*r)';
NumFinSoln = length(budget); % number of solutions returned by solver
%finSolnCount = 1;
dim = size(x0, 2); % solution dimension
if size(x0,1) >= 9
    fprintf('The number of starting solutions need to be less than 10 for RandomSearch.');
    return;
end
% Initialize
A = zeros(NumFinSoln, dim);
Afn = zeros(NumFinSoln, 1);
AFnVar = zeros(NumFinSoln, 1);

logfname = strcat(logfilename,'.txt');
logfid = fopen(logfname, 'w');

% Generate Random Variables
%structureseed = solverseed + 1;
s1 = RandStream.create('mrg32k3a', 'NumStreams', 1,'seed',solverseed);
%s2 = RandStream.create('mrg32k3a', 'NumStreams', 1,'seed',problemseed);
%s3 = RandStream.create('mrg32k3a', 'NumStreams', 1,'seed',structureseed);
RandStream.setGlobalStream(s1);


%%
% First Iteration
probStrucHandle = str2func([problem 'Structure']);
[minmax, ~, ~, ~, ~, ~, ~, ssolsV, ~, ~, ~] = probStrucHandle(numGen, problemseed);
ssolsV2 = [];
if isnan(x0)
    temp1 = randperm(numGen,s);
    ssolsV2 = ssolsV(temp1,:);
else
    temp1 = randperm(numGen,s-size(x0,1)); % Recall NumStartingSol<=9
    ssolsV2 = [x0; ssolsV(temp1,:)];
end
probHandle = str2func(problem);
[bestfn, bestfnVar, bestX] = rOptIter(probHandle, ssolsV2, s, r, problemseed, minmax);
fprintf(logfid, '======== ITERATION #%d ========\n', iterCount);
fprintf(logfid, 'x = %.4f,\n', bestX);
fprintf(logfid,'Best objective function value = %.4f,\n', bestfn);
fprintf(logfid,'Variance = %.4f.\n', bestfnVar);

for k = 2:iters(1)
    iterCount = iterCount + 1;
    problemseed = problemseed + 1;
    %call structure to figure out bounds
    [~, dim, ~, vNature, vBounds, ~, ~, ~, ~, ~, ~] = probStrucHandle(numGen, problemseed);
    %account for variable nature & generate points within boundary
    for j = 1:dim
        if vNature(j) == 0
            x = randi([vBounds(j,1), min(2^52, vBounds(j,2))]);
        else
            x = unidrnd([vBounds(j, 1), min(2^53-1, vBounds(j,2))]);
        end
        ssolsV2(j) = x;
        
    end
    
    [bestfn2, bestfnVar2, bestX2] = rOptIter(probHandle, ssolsV2, s, r, problemseed, minmax);
    if bestfn2 > bestfn
        bestfn = bestfn2;
        bestfnVar = bestfnVar2;
        bestX = bestX2;
    end
    fprintf(logfid, '======== ITERATION #%d ========\n', iterCount);
    fprintf(logfid, 'x = %.4f,\n', bestX);
    fprintf(logfid,'Best objective function value = %.4f,\n', bestfn);
    fprintf(logfid,'Variance = %.4f.\n', bestfnVar);
end
A(1,:) = bestX;
Afn(1) = bestfn;
AFnVar(1) = bestfnVar;

for nF = 2:NumFinSoln
    bestfn = -Inf;
    for k = 1:iters(nF)
        iterCount = iterCount + 1;
        problemseed = problemseed + 1;
        [~, ~, ~, VarBds, ~, ~, ~, ssolsV, ~, ~, ~] = probStrucHandle(numGen, problemseed);
        temp1 = randperm(numGen,s);
        ssolsV2 = ssolsV(temp1,:);
        [bestfn2, bestfnVar2, bestX2] = rOptIter(probHandle, ssolsV2, s, r, problemseed, minmax);
        if bestfn2 > bestfn
            bestfn = bestfn2;
            bestfnVar = bestfnVar2;
            bestX = bestX2;
        end
        fprintf(logfid, '======== ITERATION #%d ========\n', iterCount);
        fprintf(logfid, 'x = %.4f,\n', bestX);
        fprintf(logfid,'Best objective function value = %.4f,\n', bestfn);
        fprintf(logfid,'Variance = %.4f.\n\n\n', bestfnVar);
    end
    A(nF,:) = bestX;
    Afn(nF) = bestfn;
    AFnVar(nF) = bestfnVar;
end

fprintf(logfid,'Completed iterations.\n');


    function [bestfn, bestfnVar, bestX] = rOptIter(probHandle, ssolsV2, s, r, problemseed, minmax)
        fnV = zeros(s,1);  % solution objective function value & variance
        fnVarV = zeros(s,1);
        for i = 1:s
            [fn, FnVar, ~, ~, ~, ~, ~, ~] = probHandle(ssolsV2(i,:),r,problemseed);
            fnV(i) = fn;
            fnVarV(i) = FnVar;
        end
        bestfn = 0;
        bestfnInd = 0;
        if minmax == -1
            [bestfn, bestfnInd] = min(fnV);
        else
            [bestfn, bestfnInd] = max(fnV);
        end
        bestfnVar = fnVarV(bestfnInd);
        bestX = ssolsV2(bestfnInd,:);
    end

end