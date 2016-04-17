%                    SimOpt WRAPPER CODE 
%
%   ************************************************************* 
%   ***             Written by Kalyani Nagaraj                ***
%   ***         kalyanin@purdue.edu   Nov 22, 2014            ***
%   *************************************************************
%
% INPUT
%        problemname
%              Problem function name  
%        solvername
%              Solver function name
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
%
% OUTPUT
%        ncalls
%              An array (size = 'NumSoln' X 1) of budget expended 
%        x 
%              An array (size = 'NumSoln' X 'dim') of solutions
%              returned by solver
%        fn   
%              An array (size = 'NumSoln' X 1) of estimates of expected 
%              objective function value
%        FnVar
%              An array of variances/covariance matrices corresponding to  
%              the objective function at x
%              Equals NaN if solution is infeasible
%        FnGrad
%              An array of gradient estimates at x, if available
%        FnGardCov
%              An array of gradient covariance matrices at x, if available
%        constraint
%              A vector of constraint function estimators, if applicable
%        ConstraintCov
%              An array of covariance matrices corresponding to the
%              constraint function at x, if applicable
%        ConstraintGrad
%              An array of constraint gradient estimators at x, if
%              applicable
%        ConstraintGradCov
%              An array of covariance matrices of constraint gradient  
%              estimators at x, if applicable
%

function [ncalls, x, fn, FnVar, FnGrad, FnGradCov, constraint, ...
    ConstraintCov, ConstraintGrad, ConstraintGradCov]= ...
    SIMOPT(problemname, solvername, problemseed, solverseed, logfilename)

ncalls=NaN;
x=NaN;
fn=NaN;
FnVar=NaN;
FnGrad=NaN;
FnGradCov=NaN;
constraint=NaN;
ConstraintCov=NaN;
ConstraintGrad=NaN;
ConstraintGradCov=NaN;

if(numel(problemseed))~=1
    fprintf(1,'Incorrect array size of input argument "problemseed".\nSee %s.m for more information\n',problemname);
    return;
end
solverstructhandle=str2func(strcat(solvername, 'Structure'));
[NumSolverSeeds, NumStartingSols] = solverstructhandle();
if(numel(solverseed))~=NumSolverSeeds
    fprintf(1,'Incorrect array size of input argument "solverseed".\nSee %s.m for more information\n',solvername);
    return;
end
problemstructhandle=str2func(strcat(problemname, 'Structure'));
[minmax, ~, ~, ~, ~, ~, ~, x0, budget, ~, ~] ...
    = problemstructhandle(NumStartingSols, problemseed);
solverhandle=str2func(solvername);

[ncalls, x, fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ...
    ConstraintGrad, ConstraintGradCov] = ...
    solverhandle(x0, problemname, problemseed, solverseed, budget, ...
    logfilename, minmax);
end