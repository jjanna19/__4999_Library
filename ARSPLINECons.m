%==========================================================================
%                       The RSPLINE CONSTRAINED Algorithm
%==========================================================================
% DATE
%        Apr 2016
%
% AUTHOR
%        Anna Dong, Nellie Wu
% REFERENCE		
%==========================================================================
%
% INPUT
%        x0
%              Matrix (size = 'NumStartingSol' X 'dim') of 'NumStartingSol'
%              initial solutions to the solver
%              For ARSPLINECons, NumStartingSol = 1
%              Each initial solution is of size dim X 1
%        problem
%              Problem function name
%        problemseed
%              Substream index (integer >=1)
%        solverseed
%              Input seed for ARSPLINECons (integer between 1 and 2147483646)
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

%% ARSPLINECons
function [Ancalls, A, Afn, AFnVar, AFnGrad, AFnGradCov, ...
    Aconstraint, AConstraintCov, AConstraintGrad, ...
    AConstraintGradCov] = ARSPLINECons(x0, problem, problemseed, ...
    solverseed, budget, logfilename, minmax)

%% Unreported
AFnGrad = NaN;
AFnGradCov = NaN;
Aconstraint = NaN;
AConstraintCov = NaN;
AConstraintGrad = NaN;
AConstraintGradCov = NaN;




%%
% Solver parameters
kmax=1000;  % maximum number of retrospective iterations
mk=8;       % initial sample size (for k=1)
bk=10;      % mazximum number of SPLINE calls in the first retrospective 
            % iteration 
c1=1.1;     % growth rate of mk
c2=1.1;     % growth rate of bk


numfinalsoln=length(budget);
id=length(x0);

Ancalls = zeros(numfinalsoln, 1);
A = zeros(numfinalsoln, id);
Afn=zeros(numfinalsoln, 1);
AFnVar=zeros(numfinalsoln, 1);
Asiseed=zeros(numfinalsoln, 1);

% logfname=strcat(logfilename,'.txt');
% logfid=fopen(logfname, 'w');
logfid = 1;

ncalls=0;	% tracks the total calls made to the oracle
x1=x0;
k=1;
bb=1;
budg=budget(length(budget));



% Set default values.
r = 30; % Runlength time for each solution
alpha = 1; % for Nelder-Mead
gammap = 2;
betap = 0.5;
delta = 0.5; % new modification, = 0.9;
%
NumFinSoln = length(budget); % number of solutions returned by solver
dim = size(x0, 2); % solution dimension
numExtPts = dim + 1; % choonse n+1 extreme points
s1 = RandStream.create('mrg32k3a', 'NumStreams', 1,'seed',solverseed);
RandStream.setGlobalStream(s1);
problemstructhandle=str2func(strcat(problem, 'Structure'));
[minmax, ~, ~, ~, VarBds, ~, ~, xc, ~, ~, ~] = problemstructhandle(1, problemseed);
% check initSoln
if dim ~= size(xc,2)
    error('initial solution is wrong');
end
if min(budget) < r*numExtPts
    fprintf('A budget is too small for a good quality run of Nelder-Mead.');
    return
end
% calculate a scale factor ***
sf = ones(dim,2);
numL = floor(numExtPts/2);
numH = numExtPts-floor(numExtPts/2)-1;
for k12=1:dim
    if x0(k12)==VarBds(dim,1)
        numL = 0;
        numH = numExtPts-1;
    elseif x0(k12)==VarBds(dim,2)
       	numL = numExtPts;
        numH = 0;
    end
end
for k1=1:dim
    if VarBds(dim,1) > -Inf
        sf(dim,1) = -(x0(dim)-VarBds(dim,1))/numL;
    else
        sf(dim,1) = -10; %arbitrary
    end
    if VarBds(dim,2) < Inf
        sf(dim,2) = (VarBds(dim,2)-x0(dim)-10)/numH;
    else
        sf(dim,2) = 10;
    end
end
ssolsL = [x0; zeros(numL,dim)]; 
ssolsH = [x0; zeros(numH,dim)];
for k2 = 2:numL+1
    ssolsL(k2,:) = ssolsL(k2-1,:) + sf(dim,1)';
end
for k3 = 2:numH+1
    ssolsH(k3,:) = ssolsH(k3-1,:) + sf(dim,2)';
end
ssolsM = [ssolsL; ssolsH(2:end,:)]; % starting solution matrix
% used in the checkCons helper function: tstep = 0.05;

%%
iters = floor(budget/(r*numExtPts)); % stopping criterion: max # of function evals, given the problem and budget
Ancalls = (iters*r*numExtPts)';

% Initialize
A = zeros(NumFinSoln, dim);
Afn = zeros(NumFinSoln, 1);
AFnVar = zeros(NumFinSoln, 1);

logfname = strcat(logfilename,'.txt');
logfid = fopen(logfname, 'w');

probHandle = str2func(problem);

% start iterations
for nF = 1: NumFinSoln
	% Iterations
    RandStream.setGlobalStream(s1);
    problemseed2 = problemseed;
    iterCount = 1;
    
    while iterCount < iters(nF)
        fprintf(logfid, '======== ITERATION #%d ========\n', iterCount);
        [ssolsM, Plow, Flow] = iteration(ssolsM, numExtPts, r, probHandle, problemseed2, minmax, alpha, gammap, betap, delta, VarBds);
        fprintf(logfid, 'Current lowest vertex at x = %.4f,\n', Plow);
        fprintf(logfid,'Current lowest objective function value = %.4f,\n', Flow);
        iterCount = iterCount + 1;
        %problemseed2 = problemseed2 + unidrnd(100);
    end
    % Last Iteration
    [ssolsM, l2hfnVf, ~, ~, ~, l2hfnVarV] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed2, minmax, alpha,VarBds);
    A(nF,:) = ssolsM(1,:);
    Afn(nF) = -minmax*l2hfnVf(1);
    AFnVar(nF) = l2hfnVarV(1);
end

fprintf(logfid,'\nCompleted iterations.\n');

%%
% update the extreme points in one iteration
    function [ssolsMl2h, Plow, Flow] = iteration(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha, gammap, betap, delta, VarBds)
        % Iteration i
        [ssolsMl2h, l2hfnV, Prefl, Frefl, Pcent, ~] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha, VarBds);
        Plow = ssolsMl2h(1,:);
        % Accept reflection
        Flow = l2hfnV(1);
        Fsechi = l2hfnV(end-1);
        Fhigh = l2hfnV(end);
        if Frefl>=Flow && Frefl<=Fsechi
            ssolsMl2h(end,:) = Prefl;
        elseif Frefl < Flow
            Pexp2 = Prefl;
            Pexp = gammap*Prefl + (1-gammap)*Pcent;
            Pexp = checkCons(VarBds,Pexp,Pexp2);
            Fexp = -minmax*probHandle(Pexp,r,problemseed);
            if Fexp < Flow
                ssolsMl2h(end,:) = Pexp;
            else
                ssolsMl2h(end,:) = Prefl;
            end
        elseif Frefl > Fsechi
            if Frefl <= Fhigh
                ssolsMl2h(end,:) = Prefl;
                l2hfnV(end) = Frefl;
                Pcont2 = Prefl;
                Pcont = betap*Prefl + (1-betap)*Pcent;
                Pcont = checkCons(VarBds,Pcont,Pcont2);
                Fcont = -minmax*probHandle(Pcont, r, problemseed);
                if Fcont <= Fhigh
                    ssolsMl2h(end,:) = Pcont;
                else
                    for i = 2: size(ssolsMl2h,1); % from the sec extreme point to the last
                        Pnew2= Plow;
                        Pnew = delta*ssolsMl2h(i,:) + (1-delta)*Plow;
                        Pnew = checkCons(VarBds,Pnew,Pnew2);
                        ssolsMl2h(i,:) = Pnew;
                    end
                end
            end
        end
        Plow = ssolsMl2h(end,:);
        Flow = probHandle(Plow, r, problemseed);
    end



% reflect the worst point; 
% ssolsM should be constrained in the VarBds (check before inputing)
	function [ssolsMl2h, l2hfnV, Prefl, Frefl, Pcent, l2hfnVarV] = iterReflectWorst(ssolsM, numExtPts, r, probHandle, problemseed, minmax, alpha, VarBds)
        fnV = zeros(numExtPts,1);  % solution objective function value & variance
        fnVarV = zeros(numExtPts,1);
        for i = 1:numExtPts
            [fn, FnVar, ~, ~, ~, ~, ~, ~] = probHandle(ssolsM(i,:),r,problemseed); 
            fnV(i) = -minmax*fn; % convert maximization problem to minimization -objective fcn problem
            fnVarV(i) = FnVar;
        end
        [l2hfnV,l2hfnIndV] = sort(fnV);
        l2hfnVarV = fnVarV(l2hfnIndV,:);
        
        Phigh = ssolsM(l2hfnIndV(end),:);

        ssolsMl2h = ssolsM(l2hfnIndV,:); 
        Pcent = mean(ssolsMl2h(1:end-1,:));
        Prefl2 = Phigh; %save the original point
        Prefl = (1+alpha)*Pcent - alpha*Phigh;
        Prefl = checkCons(VarBds,Prefl,Prefl2); %check if the new reflect point respects VarBds (if not, change it)
        
        [Frefl, FreflVar, ~, ~, ~, ~, ~, ~] = probHandle(Prefl,r,problemseed);
        Frefl = -minmax*Frefl;
    end


% check and modify (if needed) the new point, based on VarBds.
    function modiSsolsV = checkCons(VarBds,ssolsV,ssolsV2) 
        tstep = 0.05; %arbitrary***
        col = size(ssolsV,2);
        stepV = ssolsV - ssolsV2;
        t = 1; % t>0 for the correct direction
        for j = 1:col
            while ssolsV(j)<VarBds(j,1) || ssolsV(j)>VarBds(j,2)
                t = t-tstep;
                ssolsV = ssolsV2+t*stepV;
            end
        end
        modiSsolsV = ssolsV;
    end


end

