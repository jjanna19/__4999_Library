%==========================================================================
%                        The R-SPLINE Algorithm 
%==========================================================================
% DATE
%        November 2014
%
% AUTHOR
%        Kalyani Nagaraj
%        kalyanin AT purdue DOT edu
%
% REFERENCE
%        H. Wang, B. Schmeiser, and R. Pasupathy. ACM TOMACS. 2013
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
%              Input seed for R-SPLINE (integer between 1 and 2147483646)
%              See function u16807d for details.
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

%% RSPLINE
function [Ancalls, A, Afn, AFnVar, AFnGrad, AFnGradCov, ...
     Aconstraint, AConstraintCov, AConstraintGrad, ...
     AConstraintGradCov] = RSPLINE(x0, problem, problemseed, ...
     solverseed, budget, logfilename, minmax)
 
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
AFnGrad=NaN;
AFnGradCov=NaN;
Aconstraint=NaN;
AConstraintCov=NaN;
AConstraintGrad=NaN;
AConstraintGradCov=NaN;
Asiseed=zeros(numfinalsoln, 1);

% logfname=strcat(logfilename,'.txt');
% logfid=fopen(logfname, 'w');
logfid = 1;

ncalls=0;	% tracks the total calls made to the oracle
x1=x0;
k=1;
bb=1;
budg=budget(length(budget));

%Begin Retrospective Iterations
while ncalls<budg && k<=kmax

    % fprintf(logfid, '==== RETROSPECTIVE ITERATION k = %d ====\n', k);
    % For deterministic constraints
    mk=ceil(mk*c1);
    bk=ceil(bk*c2);
    
    % fprintf(logfid, 'mk = %d, bk = %d\n', mk, bk);
    
    xk=x1;

    % fprintf(logfid, '===== BEGIN SPLINE =====\n');

	[splinencalls, x1, x1fn, x1FnVar, ~, ~, ...
        solverseed] = SPLINE(problem, xk, mk, bk, problemseed, solverseed, ...
        logfid, minmax);	       
    if splinencalls == -1     % initial solution is infeasible!
		return
    end
    ncalls = ncalls + splinencalls; %splinencalls is the number of oracle 
                                    %calls for each call of SPLINE

    if ncalls<=budget(bb)
        Ancalls(bb)=ncalls;
        A(bb,:) = x1';
        Afn(bb)=-minmax*x1fn; %Afn(bb)=x1fn;
        AFnVar(bb)=x1FnVar;
        Asiseed(bb)=solverseed;
    elseif bb<length(budget)
        bb=bb+1;
        Ancalls(bb)=ncalls;
        A(bb,:) = x1';
        Afn(bb)=-minmax*x1fn; %Afn(bb)=x1fn;
        AFnVar(bb)=x1FnVar;
        Asiseed(bb)=solverseed;
    end
    
    
        % fprintf(logfid, '===== SPLINE ENDED =====\nncalls so far = %d, xbest.ghat = %0 .12f, xbest = [', ncalls, x1fn);
        for i=1:id
			% fprintf(logfid, '%d ', x1(i));	
        end
        % fprintf(logfid, ']\n');
        
        
	k=k+1;
end
% fclose(logfid);
end

%% SPLINE
function [ncalls, xnew, xnewfn, xnewFnVar, xnewconstraint, ...
    xnewConstraintCov, sseed]=SPLINE(problem, x0, mk, bk, iseed, sseed, ...
    logfid, minmax)
		
id=length(x0);
ncalls=0;
xnew=x0;

[xnewfn, xnewFnVar, ~, ~, xnewconstraint, xnewConstraintCov, ~, ~] = ...
    oracle(str2func(problem), xnew, mk, iseed, minmax);

if isnan(xnewfn) % infeasible solution
    ncalls=-1;
    return
end
x0fn=xnewfn;
x0FnVar=xnewFnVar;
x0constraint=xnewconstraint;
x0ConstraintCov=xnewConstraintCov;

% fprintf(logfid, '\t===== BEGIN SPLINE LOOP ===\n');
for i=1:bk
	% fprintf(logfid, '\t\t=== bk = %d ===\n', i);
		
	[SPLIncalls, xold, xoldfn, xoldFnVar, xoldcontsraint, ...
        xoldConstraintCov, sseed] = SPLI(problem, xnew, xnewfn, ...
        xnewFnVar, xnewconstraint, xnewConstraintCov, mk, iseed, sseed, ...
        logfid, minmax);
	
		% fprintf(logfid, '\t\t\t==== SPLI ends\n\t\t\txold = [');	
        for j=1:id
			% fprintf(logfid, '%d ', xold(j));
        end
		% fprintf(logfid, ']\n\t\t\txold.ghat = %.12f\n', xoldfn);
		
	[NEncalls, xnew, xnewfn, xnewFnVar, xnewconstraint, ...
        xnewConstraintCov] = NE(problem, xold, xoldfn, ...
        xoldFnVar, xoldcontsraint, xoldConstraintCov, mk, iseed, ...
        logfid, minmax);

        % fprintf(logfid, '\t\t\t==== NE ends\n\t\t\txnew = [');	
        for j=1:id
            % fprintf(logfid, '%d ', xnew(j));
        end
        % fprintf(logfid, ']\n\t\t\txnew.ghat = %.12f\n', xnewfn);

    ncalls = ncalls + SPLIncalls + NEncalls;

        % fprintf(logfid, '\n\n\t\tSPLINE ncalls = prev ncalls + SPLI calls + NE calls = %d\n', ncalls);
        
    if  xoldfn==xnewfn
		% fprintf(logfid, '\n\t\tSPLINE ended at bk=%d since NE and SPLI returned the same solution\n\n', i);
		break
    end
end    
	
% starting solution is better than the new solution by a very small margin
if x0fn + 0.00005 <= xnewfn
	xnew=x0;
	xnewfn=x0fn;
    xnewFnVar=x0FnVar;
    xnewconstraint=x0constraint;
    xnewConstraintCov=x0ConstraintCov;
end
end

%% NE
function [ncalls, xnew, xnewfn, xnewFnVar, xnewconstraint, ...
    xnewConstraintCov] = NE(problem, x, fn, FnVar, constraint, ...
    ConstraintCov, mk, iseed, logfid, minmax)

id = length(x); 
xold=x;
xnew=x;
xnewfn=fn;
xnewFnVar=FnVar;
xnewconstraint=constraint;
xnewConstraintCov=ConstraintCov;

ncalls=0;
y2=fn;
ixquad=zeros(id,1);

	% fprintf(logfid, '\n\t\t\t===== NE BEGINS =====\n');
	% fprintf(logfid, '\t\t\tghat at center = %.12f, center = [', y2);
    for i=1:id
		% fprintf(logfid, '%d ', xnew(i));
    end
	% fprintf(logfid, ']\n');

for i=1:id
	count=1;

    xold(i)=xold(i)+1;
    [xoldfn, xoldFnVar, ~, ~, xoldconstraint, xoldConstraintCov, ~, ~]...
        =oracle(str2func(problem), xold, mk, iseed, minmax);

    if ~isnan(xoldfn)
		ncalls = ncalls + mk;
		y1=xoldfn; 
		count=count+1;
        if xoldfn + 0.00005 < xnewfn
			xnew=xold;
			xnewfn=xoldfn;
			xnewFnVar=xoldFnVar;
            xnewconstraint=xoldconstraint;
            xnewConstraintCov=xoldConstraintCov;
        end
    end

	xold(i)=xold(i)-2;
    [xoldfn, xoldFnVar, ~, ~, xoldconstraint, xoldConstraintCov, ~, ~]=...
        oracle(str2func(problem), xold, mk, iseed, minmax);
    if ~isnan(xoldfn)
		ncalls = ncalls + mk;
		y3=xoldfn; 
		count=count+1;
        if xoldfn + 0.00005 < xnewfn
			xnew=xold;
			xnewfn=xoldfn;
            xnewFnVar=xoldFnVar;
            xnewconstraint=xoldconstraint;
            xnewConstraintCov=xoldConstraintCov;
        end 
    end
	xold(i)=xold(i)+1;
	xqnew=xold(i);

    %quadratic search
    if count==3 
		a = (y1+y3)/2.0 - y2;
		b = (y1-y3)/2.0;
        if a-0.00005 > 0
			xqnew = int32(xold(i) - (b / (a + a)));
        end
		% fprintf(logfid, '\t\t\t\ti = %d, a = %.12f, b = %.12f, xqnew = %.12f\n', i, a, b, xqnew);
		% fprintf(logfid, '\t\t\t\ty2 = %.12f, y1 = %.12f, y3 = %.12f\n', y2, y1, y3);
    end
    if  abs(xqnew) < 2147483646.0 %2^31-2
		ixquad(i) = xqnew;
    end
	% fprintf(logfid, '\t\t\t\txold[%d] = %d, ixquad[%d] = %d\n', i, xold(i), i, ixquad(i));
end
		
[ixquadfn, ixquadFnVar, ~, ~, ixquadconstraint, ixquadConstraintCov, ...
    ~, ~] = oracle(str2func(problem), ixquad, mk, iseed, minmax);

if ~isnan(ixquadfn)
	ncalls = ncalls + mk; 
    if ixquadfn + 0.00005 < xnewfn
		xnew=ixquad;
        xnewfn=ixquadfn;
        xnewFnVar=ixquadFnVar;
        xnewconstraint=ixquadconstraint;
        xnewConstraintCov=ixquadConstraintCov;
    end
end

        % fprintf(logfid, '\t\t\tixquad.ghat = %.12f, ixquad = [', ixquadfn);
        for i=1:id
			% fprintf(logfid, '%d ', ixquad(i));
        end
		% fprintf(logfid, ']\n');
end


%% SPLI
function [ncalls, xbest, xbestfn, xbestFnVar, xbestconstraint, ...
    xbestConstraintCov, sseed] = SPLI(problem, x, fn, FnVar, ...
    constraint, ConstraintCov, mk, iseed, sseed, logfid, minmax)

id = length(x); 
imax=100;
jmax=5;

% fprintf(logfid, '\t\t\tBEGIN SPLI ===\n');

xbest=x;
xbestfn=fn;
xbestFnVar=FnVar;
xbestconstraint=constraint;
xbestConstraintCov=ConstraintCov;

ncalls=0;
s0=2.0;
c=2.0;
			
for j=0:jmax
    %PRINT TO LOGFILE
    % fprintf(logfid, '\t\t\t\t=== j = %d ===\n', j);
    % fprintf(logfid, '\t\t\t\ti/p seed to PERTURB = %d\n', sseed);
    %END PRINT

    [x1, sseed]=PERTURB(xbest, sseed, logfid);
			
    %PRINT TO LOGFILE
    % fprintf(logfid, '\t\t\t\tPertrubed value of xbest = [');
    for i=1:id
        % fprintf(logfid, '%0.12f ', x1(i));
    end
    % fprintf(logfid, ']\n');
    % fprintf(logfid, '\t\t\t\tPLI begins ===\n');		
    %END PRINT

   [fbar, gamma, npoints, plixbest, plixbestfn, plixbestFnVar, ...
       plixbestconstraint, plixbestConstraintCov] = PLI(problem, ...
       x1, mk, iseed, logfid, minmax);
		
            % fprintf(logfid, '\t\t\t\t=== PLI ends\n');
			% fprintf(logfid, '\t\t\t\tgbar = %.12f\n\t\t\t\tgamma = [',fbar);
            for i=1:id
				% fprintf(logfid, '%.12f ', gamma(i));
            end
			% fprintf(logfid, ']\n\t\t\t\txbest.ghat = %.12f\n\t\t\t\txbest = [',plixbestfn);
            for i=1:id	
				% fprintf(logfid, '%d ', plixbest(i));
            end
            % fprintf(logfid,']\n\t\t\t\tnpoints=%d\n', npoints);

	ncalls = ncalls + npoints*mk;	

	%regardless of whether npoints=id+1 or not, update current best
    if  plixbestfn + 0.00005 < xbestfn && npoints>0
        xbest=plixbest;
        xbestfn=plixbestfn;
        xbestFnVar=plixbestFnVar;
        xbestconstraint=plixbestconstraint;
        xbestConstraintCov=plixbestConstraintCov;
    end
		
    if npoints < id+1
        return 
    end
			
	glength=norm(gamma);
	
	% fprintf(logfid, '\t\t\t\tglength = %.12f\n', glength);
		
    if glength + 0.00005 <= 0
		return
    end
    x0=xbest;
    gamma=gamma/glength;

    for i=0:imax
        s = s0 * c^i;
        ix1=zeros(id,1);
        for k=1:id
			ix1(k)=floor(x0(k)-s*gamma(k)+0.5);
        end
        [ix1fn, ix1FnVar, ~, ~, ix1constraint, ix1ConstraintCov, ~, ~] ...
            = oracle(str2func(problem), ix1, mk, iseed, minmax);

        % PRINT TO LOGFILE
        % fprintf(logfid, '\t\t\t\t\t== i = %d ==\n\t\t\t\t\ts = %.12f, ix1 = [', i, s);
        for k=1:id 
            % fprintf(logfid, '%d ', ix1(k));
        end
        % fprintf(logfid, ']\n\t\t\t\t\tix1.ghat = %.12f\n', ix1fn);
        % END PRINT

        if isnan(ix1fn) % if ix1 is infeasible
            return	
        end
		ncalls = ncalls + mk;
        if ix1fn >= xbestfn + 0.00005 && i <= 2
            return
        end
        if ix1fn >= xbestfn + 0.00005
            break
        end
        
		xbest=ix1;
		xbestfn = ix1fn;
		xbestFnVar=ix1FnVar;
        xbestconstraint=ix1constraint;
        xbestConstraintCov=ix1ConstraintCov;
    end    
end
end


%% PLI
function [fbar, gamma, npoints, plixbest, plixbestfn, plixbestFnVar, ...
    plixbestconstraint, plixbestConstraintCov] = PLI(problem, x, ...
    mk, iseed, logfid, minmax)

id = length(x); npoints=0;
strange=3.145962987654;
gamma=zeros(id,1);
x0=floor(x);
%z=zeros(id,1);
z=x-x0;
z=[1;z;0];

[~, p]=sort(z, 'descend');
w=zeros(id+1,1);
for i=1:id+1		
	w(i,1)=z(p(i))-z(p(i+1));
end	
wsum=0;
fbar=0;

%PRINT TO LOGFILE
% fprintf(logfid, '\t\t\t\t\tp z w\n');
for i=1:id+1
    % fprintf(logfid, '\t\t\t\t\t%d %0.15f %0.15f\n', p(i), z(i), w(i));
end
% fprintf(logfid, '\n\t\t\t\t\ti=1, mk=%d\n', mk);	
% fprintf(logfid, '\t\t\t\t\tx0 = [');
for i=1:id
    % fprintf(logfid, '%d ', x0(i));
end
% fprintf(logfid, ']\n');
%END PRINT

[x0fn, x0FnVar, ~, ~, x0constraint, x0ConstraintCov, ~, ~]=...
    oracle(str2func(problem), x0, mk, iseed, minmax);	
if ~isnan(x0fn) 
	npoints=npoints+1;
	wsum = wsum + w(1);
	fbar = fbar + w(1)*x0fn;
	ghatold = x0fn;
    plixbest=x0;
    plixbestfn = x0fn;
    plixbestFnVar=x0FnVar;
	plixbestconstraint=x0constraint;
	plixbestConstraintCov=x0ConstraintCov;	
    
    % fprintf(logfid, '\t\t\t\t\tx0.ghat = %.12f, w = %.12f, ghatold = %.12f, gbar = %.6f\n\n', x0fn, w(1), ghatold, fbar);
else
    ghatold = 0;
	plixbestfn = strange;
end

for i=2:id+1
    x0(p(i)-1)=x0(p(i)-1)+1;
    
    % PRINT TO LOGFILE
    % fprintf(logfid, '\n\t\t\t\t\ti = %d', i);
    % fprintf(logfid, ', mk = %d\n', mk);		
    % fprintf(logfid, '\t\t\t\t\tx0 = [');
    for j=1:id
        % fprintf(logfid, '%d ', x0(j));
    end
    % fprintf(logfid, ']\n');
    % END PRINT
    
    %call oracle at the other id points that form the simplex
	[x0fn, x0FnVar, ~, ~, x0constraint, x0ConstraintCov, ~, ~]=...
        oracle(str2func(problem), x0, mk, iseed, minmax);
    if ~isnan(x0fn)
		npoints=npoints+1;
		wsum = wsum + w(i);
		fbar = fbar + w(i)*x0fn;
		gamma(p(i)-1) = x0fn - ghatold;
		ghatold = x0fn;

        %PRINT TO LOGFILE
        % fprintf(logfid, '\t\t\t\t\tx0.ghat = %.12f, w = %.12f, ghatold = %.12f\n', x0fn, w(i), ghatold);
        % fprintf(logfid, '\t\t\t\t\tgbar=%0.6f, npoints=%d\n', fbar, npoints);
        %END PRINT		
        
        if plixbestfn == strange || x0fn + 0.00005 < plixbestfn
            plixbest=x0;
            plixbestfn = x0fn;
            plixbestFnVar=x0FnVar;
            plixbestconstraint=x0constraint;
            plixbestConstraintCov=x0ConstraintCov;	
        end
        
        %PRINT TO LOGFILE
        % fprintf(logfid, '\t\t\t\t\txbest.ghat = %.12f\n\t\t\t\t\txbest = [', plixbestfn);
        for j=1:id
            % fprintf(logfid, '%d ', plixbest(j));
        end
        % fprintf(logfid, ']\n');
        %END PRINT
    end
end

if wsum > 0.00005
    fbar = fbar/wsum;	
end

%PRINT TO LOGFILE
% fprintf(logfid, '\n\t\t\t\t\tgbar = %.12f, npoints = %d\n', fbar, npoints);
%END PRINT

end


%% PERTURB
function [xpert, sseed] = PERTURB(x, sseed, logfid)

id=length(x);
xpert=zeros(id,1);
	% fprintf(logfid, '\t\t\t\tPERTURB begins ==== \n');
for i=1:id
        % fprintf(logfid, '\t\t\t\t\t sseed = %d, ', sseed);
	[sseed, u] = u16807d(sseed);
    xpert(i) = x(i) + .3*(u - 0.5);
    % fprintf(logfid, 'u = %0.12f, xpert[%d] = %.12f\n', u, i, xpert(i));
end
    % fprintf(logfid, '\t\t\t\t==== PERTURB ends\n');
end


%% u16807d
function [iseed,u16807d]=u16807d(iseed)
%..........................................................................
%     bruce schmeiser     january 1989.                                   .
%     a linear congruential pseudorandom number generator                 .
%       using constant 16807 and modulus (2**31)-1.                       .
%              iseed = iseed*16807 (mod 2^31 -1)                          .
%     for implementations that don't require double precision, see        .
%       s.k. park and k.w. miller, "random numbers generators: good       .
%         ones are hard to find," cacm, 31, 10 (October 1988), 1192-1201. .
%     in correct implementations, starting with a seed value of 1 will    .
%       result in a seed value of 1043618065 on call number 10001.        .
%..........................................................................
%     input:  iseed.   integer.                                           . 
%                        chosen from [1,2147483646] on the first call.    .
%                        thereafter, the value returned from the last call.
%     output: iseed.   integer.                                           .
%                        to be used in the next call.                     .
%     output: u16807d.  real.                                             .
%                        a pseudorandom number in (0,1).                  .
%..........................................................................
u16807d=0;
while (u16807d<=0 || u16807d>=1)
    iseed = mod (iseed * 16807,2147483647);
    u16807d = iseed / 2147483648;
end
end


%% oracle
function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ...
    ConstraintGrad, ConstraintGradCov] = oracle(orchandle, x, m, iseed, minmax)

%
% Input: 
%  orchandle
%    problem function handle
%  x
%    decision variable (size = 1 X id)
%  m 
%    sample size
%  iseed
%    input seed substream index (integer >=1)
% 
% Output:
%  fn
%    estimate of expected real valued objective function 
%    equals NaN if solution is infeasible
%  FnVar
%    variance of expected objective function
%    equals NaN if solution is infeasible
%  FnGrad
%    vector of estimate of gradient at x
%  FnGardCov
%    covariance matrix of gradient estimator at x
%  constraint
%    vector of constraint function estimators
%  ConstraintCov
%    covariance matrix of constraint extimators at x
%  ConstraintGrad
%    matrix of gradient estimators of constraint functions at x
%  ConstraintGradCov
%    3-dimensional covariance array of gradient estimators of constraint
%    functions at x

[fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ...
    ConstraintGrad, ConstraintGradCov] = orchandle(x, m, iseed, NaN);
fn=-minmax*fn;
end