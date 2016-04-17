%==========================================================================
%                       The Nelder-Mead Algorithm
%==========================================================================
% DATE
%        Feb 2016
%
% AUTHOR
%        Anna Dong, Nellie Wu
% REFERENCE		Russell R. Barton, John S. Ivey, Jr., (1996)
%		  			Nelder-Mead Simplex Modifications for Simulation 
%			 		Optimization. Management Science 42(7):954-973.
%==========================================================================
%
% INPUT
%        x0
%              Matrix (size = 'NumStartingSol' X 'dim') of 'NumStartingSol'
%              initial solutions to the solver
%              For Nelder-Mead, NumStartingSol<=9
%              Each initial solution is of size dim X 1
%        problem
%              Problem function name
%        problemseed
%              Substream index (integer >=1)
%        solverseed
%              Input seed for Nelder-Mead (integer between 1 and 2147483646)
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

%% Nelder-Mead
function [Ancalls, A, Afn, AFnVar, AFnGrad, AFnGradCov, ...
    Aconstraint, AConstraintCov, AConstraintGrad, ...
    AConstraintGradCov] = NM(x0, problem, problemseed, ...
    solverseed, budget, logfilename, minmax)
%% Unreported
AFnGrad = NaN;
AFnGradCov = NaN;
Aconstraint = NaN;
AConstraintCov = NaN;
AConstraintGrad = NaN;
AConstraintGradCov = NaN;

% Set default values.
r = 30; % Runlength time for each solution
extP = 0.05; % Generate extreme points around x0 by adding 5% of each x0(k) to x0
ext0 = 0.00025; % use ext0 for 0 values
alpha = 1; % for Nelder-Mead
gammap = 2;
betap = 0.5;
delta = 0.5; % new modification, = 0.9;

%%
NumFinSoln = length(budget); % number of solutions returned by solver
dim = size(x0, 2); % solution dimension
numExtPts = dim + 1; % choonse n+1 extreme points
% ssolsM = [x0; zeros(dim)]; % starting solution matrix
% for i = 1:dim
%     temp1 = x0;
%     temp1(i) = x0(i)*(1+extP);
%     if x0(i)==0
%         temp1(i) = ext0;
%     end
%     ssolsM(i+1,:) = temp1;
% end
ssolsM = [x0; zeros(dim)]; % starting solution matrix
for i = 1:dim
	for k = 1:dim
		ssolsM(i+1,k) = ssolsM(i,k)*(1+extP);
		if ssolsM(i,k) == 0
			ssolsM(i+1,k) = ext0;
		end
	end
end

if min(budget) < r*numExtPts
    fprintf('A budget is too small for a good quality run of Nelder-Mead.');
    return
end
iters = floor(budget/(r*numExtPts)); % stopping criterion: max # of function evals, given the problem and budget
iterCount = 1; % which budget we are currently looking at
Ancalls = (iters*r*numExtPts)';

% Initialize
A = zeros(NumFinSoln, dim);
Afn = zeros(NumFinSoln, 1);
AFnVar = zeros(NumFinSoln, 1);

logfname = strcat(logfilename,'.txt');
logfid = fopen(logfname, 'w');

probHandle = str2func(problem);
for nF = 1: length(iters)
	% Iterations
	while iterCount < iters(nF)
		fprintf(logfid, '======== ITERATION #%d ========\n', iterCount);
		[ssolsM, Plow, Flow] = iteration(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha, gammap, betap, delta);
		fprintf(logfid, 'Current lowest vertex at x = %.4f,\n', Plow);
		fprintf(logfid,'Current lowest objective function value = %.4f,\n', Flow);
		iterCount = iterCount + 1;
		problemseed = 0.5*problemseed + 0.5*solverseed + 1; % or problemseed = problemseed + 1, no use of solverseed
	end
	% Last Iteration
	[ssolsM, l2hfnVf, ~, ~, ~, l2hfnVarV] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha);
	A(nF,:) = ssolsM(1,:);
    Afn(nF) = l2hfnVf(1);
    AFnVar(nF) = l2hfnVarV(1);
    % reset iterCount
	iterCount = 1;
    if minmax == 1
	Afn(nF) = -Afn(nF);
    end
end

fprintf(logfid,'\nCompleted iterations.\n');



function [ssolsMl2h, Plow, Flow] = iteration(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha, gammap, betap, delta)
	% Iteration i
	[ssolsMl2h, l2hfnV, Prefl, Frefl, Pcent, ~] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha);
	Plow = ssolsMl2h(1,:);
	% Accept reflection
	Flow = l2hfnV(1);
	Fsechi = l2hfnV(end-1);
    Fhigh = l2hfnV(end);
	if Frefl>=Flow && Frefl<=Fsechi
		ssolsMl2h(end,:) = Prefl;
	elseif Frefl < Flow
	    Pexp = gammap*Prefl + (1-gammap)*Pcent;
	    Fexp = probHandle(Pexp,r,problemseed);
	    if Fexp < Flow 
	        ssolsMl2h(end,:) = Pexp;
	    else 
	    	ssolsMl2h(end,:) = Prefl;
	    end
	elseif Frefl > Fsechi
		if Frefl <= Fhigh 
			ssolsMl2h(end,:) = Prefl;
			l2hfnV(end) = Frefl;
			Pcont = betap*ssolsMl2h(end,:) + (1-betap)*Pcent;
			[Fcont, ~, ~, ~, ~, ~, ~, ~] = probHandle(Pcont, r, problemseed);
	    	if Fcont <= Fhigh
	        	ssolsMl2h(end,:) = Pcont; 
	    	else
				for i = 2: size(ssolsMl2h,1); % from the sec extreme point to the last
	                ssolsMl2h(i,:) = delta*ssolsMl2h(i,:) + (1-delta)*Plow;
	            end
	        end
	    end	
	end
end



	function [ssolsMl2h, l2hfnV, Prefl, Frefl, Pcent, l2hfnVarV] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha)
        fnV = zeros(numExtPts,1);  % solution objective function value & variance
        fnVarV = zeros(numExtPts,1);
        for i = 1:numExtPts
            [fn, FnVar, ~, ~, ~, ~, ~, ~] = probHandle(ssolsM(i,:),r,problemseed);
            fnV(i) = fn;
            fnVarV(i) = FnVar;
        end
        if minmax == 1	% convert maximization problem to minimization -objective fcn problem
        	fnV = -fnV;
        	fnVarV(i) = FnVar;
        end
        [l2hfnV,l2hfnIndV] = sort(fnV);
        l2hfnVarV = fnVarV(l2hfnIndV,:);

        Phigh = ssolsM(l2hfnIndV(end),:);
        Psechi = ssolsM(l2hfnIndV(end-1),:);
        Plow = ssolsM(l2hfnIndV(1),:);
        Fhigh = fnV(end);
        Fsechi = fnV(end-1);
        Flow = fnV(1);

        ssolsMl2h = ssolsM(l2hfnIndV,:);
        Pcent = sum(ssolsMl2h(1:end-1,:));
        Prefl = (1+alpha)*Pcent - alpha*Phigh;
        [Frefl, FreflVar, ~, ~, ~, ~, ~, ~] = probHandle(Prefl,r,problemseed);
	end

end