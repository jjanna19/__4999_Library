%==========================================================================
%                The More Conservative Random Search Algorithm
%==========================================================================
% DATE
%        Feb 2016
%
% AUTHOR
%        Anna Dong
%
%==========================================================================
%
% INPUT
%        x0
%              Matrix (size = 'NumStartingSol' X 'dim') of 'NumStartingSol'
%              initial solutions to the solver, can be NaN
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
    AConstraintGradCov] = RSCons(x0, problem, problemseed, ...
    solverseed, budget, logfilename, minmax)
%tic
%% Unreported
AFnGrad = NaN;
AFnGradCov = NaN;
Aconstraint = NaN;
AConstraintCov = NaN;
AConstraintGrad = NaN;
AConstraintGradCov = NaN;

%%
r = 30; % Runlength time for each solution
if min(budget) < size(x0,1)*r
    fprintf('A budget is too small to run all starting solutions.');
    return;
end
NumFinSoln = length(budget); % number of solutions returned by solver
numAttemptsV = floor(budget/r); 
numGenV = numAttemptsV*100; % Size of the pool to (randomly) search among at each iteration
dim = size(x0, 2); % solution dimension

Ancalls = (numAttemptsV*r)';

% Initialize
A = zeros(NumFinSoln, dim);
Afn = zeros(NumFinSoln, 1);
AFnVar = zeros(NumFinSoln, 1);

logfname = strcat(logfilename,'.txt');
logfid = fopen(logfname, 'w');

% Generate Random Variables
s1 = RandStream.create('mrg32k3a', 'NumStreams', 1,'seed',solverseed);


%%
probStrucHandle = str2func([problem 'Structure']);
probHandle = str2func(problem);
for k = 1:NumFinSoln
    RandStream.setGlobalStream(s1);
    [minmax, ~, ~, ~, ~, ~, ~, ssolsV, ~, ~, ~] = probStrucHandle(numGenV(k), problemseed);
    ssolsV2 = [];
    if isnan(x0)
        temp1 = randperm(numGenV(k),numAttemptsV(k));
        ssolsV2 = ssolsV(temp1,:);
    else
        temp1 = randperm(numGenV(k),numAttemptsV(k)-size(x0,1)); % run starting solns
        ssolsV2 = [x0; ssolsV(temp1,:)];
    end
    [bestfn, bestfnVar, bestX] = findOpt(probHandle, ssolsV2, r, problemseed, minmax);
    A(k,:) = bestX;
    Afn(k) = bestfn;
    AFnVar(k) = bestfnVar;
end
% fprintf(logfid, '======== Budget #%d ========\n', k);
% fprintf(logfid, 'x = %.4f,\n', bestX);
% fprintf(logfid,'Best objective function value = %.4f,\n', bestfn);
% fprintf(logfid,'Variance = %.4f.\n', bestfnVar);


    function [bestfn, bestfnVar, bestX] = findOpt(probHandle, ssolsV2, r, problemseed, minmax)
        s = length(ssolsV2);
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
%toc
end