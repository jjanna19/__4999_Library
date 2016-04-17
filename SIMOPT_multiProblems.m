%                    SimOpt WRAPPER CODE MULTIPROBLEMS
%
%   ************************************************************* 
%   ***             Written by Anna Dong, Nellie Wu           ***
%   ***                       Mar 12, 2016                    ***
%   *************************************************************
%
% INPUT
%        problemArray
%              Problem function numbers (as an array)  
%        solvernumber
%              Solver function number
%        problemseed
%              Index of the substreams to use (integer >= 1)
%              See problem code for more information.
%        solverseed
%              Input seed(s) for solver, if required. 
%              Should equal NaN otherwise. 
%        logfilename
%              Character string indicating path to a log file (if
%              required)
%              Should equal NaN otherwise
%        budgetV: vector of budgets
%        plotopt: option to plot or not, default = 1, yes
%
% OUTPUT
%        Plot: x-axis: budget, 
%               y-axis: objective value, legend for each problem
%        FMatrix: length(problemArray) X length(budget)
%

function FMatrix = SIMOPT_multiProblems(problemArray, solvernumber, problemseed, solverseed, budgetV, logfilename, plotopt)
%initialize
FMatrix=zeros(length(problemArray),length(budgetV));
FVarMatrix=zeros(length(problemArray),length(budgetV));
close all
hold on


%input
repsAlg = 30;
reps = 30;
CILevel = 0.95;
ProblemStrV0 = char(importdata(strcat(pwd, '/index/ProblemNumber.csv')));
SolverStrV0 = char(importdata(strcat(pwd, '/index/SolverNumber.csv')));

%
problemnameV = repmat('', length(problemArray), size(ProblemStrV0,2));
solvername = strtrim(SolverStrV0(solvernumber,:));
solverstructhandle=str2func(strcat(solvername, 'Structure'));
[NumSolverSeeds, NumStartingSols] = solverstructhandle();
solverhandle=str2func(solvername);

    for i=1:length(problemArray)
        FMatrix2 = zeros(repsAlg,length(budgetV));
        problemnameV(i,:) = ProblemStrV0(problemArray(i),:);
        problemname = strtrim(problemnameV(i,:));
        problemstructhandle=str2func(strcat(problemname, 'Structure'));
        [minmax, ~, ~, ~, ~, ~, ~, x0, ~, ~, ~] ...
        = problemstructhandle(NumStartingSols, problemseed);
        problemhandle=str2func(problemname);
        for j=1:repsAlg  %do 30 replications for each algo on the problem
            [minmax, ~, ~, ~, ~, ~, ~, x0, ~, ~, ~] = problemstructhandle(NumStartingSols, problemseed);
            [~, soln, ~, ~, ~, ~, ~, ~, ~, ~] = ...
            solverhandle(x0, problemname, problemseed, solverseed, budgetV, logfilename, minmax);
            seed1=unidrnd(10000,1,1);           
            for k=1:size(soln,1)
                [fn, ~, ~, ~, ~, ~, ~, ~] = problemhandle(soln(k,:), reps, seed1);
                FMatrix2(j,k) = fn;
            end
        end
        FMatrix(i,:) = mean(FMatrix2);
        FVarMatrix(i,:) = var(FMatrix2);
        FnSEM = sqrt(FVarMatrix(i,:)/repsAlg);
        EWidth = norminv(1-(1-CILevel)/2,0,1)*FnSEM;
        if plotopt == 1
            rng('shuffle');
            errorbar(budgetV,FMatrix(i,:),EWidth,'-o', 'Color',rand(1,3));
        end
    end
    if plotopt == 1
        legend(problemnameV);
        xlabel('Budget'); ylabel('Objective value');
        title(['Solver: ',solvername]);
    end
    
end
