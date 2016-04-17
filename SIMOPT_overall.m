%                    SimOpt WRAPPER CODE OVERALL
%
%   ******************************************************************* 
%   *** By Anna Dong, Nellie Wu, Arielle Anderer, Jennifer Mallette ***
%   ***                         Mar 12, 2016                        ***
%   *******************************************************************
%
% INPUT
%        problemArray
%              Problem function number   
%        solverArray
%              Solver function numbers (as an array)
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
%        budget: vector of budgets
%
% OUTPUT
%        FMatrix: length(problemArray) X length(budget)
%

function FMatrix = SIMOPT_overall(problemArray, solverArray, problemseed, solverseed, budgetV, logfilename)
%initialize
FMatrix=zeros(length(problemArray),length(solverArray)*length(budgetV));

%input
ProblemStrV0 = char(importdata(strcat(pwd, '/index/ProblemNumber.csv')));
SolverStrV0 = char(importdata(strcat(pwd, '/index/SolverNumber.csv')));

%
for i=1:length(problemArray)
    problemnumber = problemArray(i);
    FMatrix2 = SIMOPT_multiSolvers(problemnumber, solverArray, problemseed, solverseed, budgetV, logfilename, 0);
    FMatrix(i,:) = reshape(FMatrix2',[1,length(solverArray)*length(budgetV)]);
end

save(strcat(pwd, '/', '/output/Evaluation_Overall.csv'), 'FMatrix', '-ascii');

end