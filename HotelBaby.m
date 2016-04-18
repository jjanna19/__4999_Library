%   ***************************************
%   *** Code written by German Gutierrez***
%   ***         gg92@cornell.edu        ***
%   ***************************************

% RETURNS: Total revenue captured
function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = HotelBaby(x, runlength, seed, ~)
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = HotelBaby(x, runlength, seed, other);
% x is a vector of reservations for each period, in order [r2 r3...rT 0];
% runlength is the number of days of demand to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used

FnGrad = NaN;
FnGradCov = NaN;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;

C=100;
if (sum(x < 0)>0) || sum(x>C)>0  || (runlength <= 0) || (runlength ~= round(runlength)) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('x should be >= 0 and <= C, runlength should be positive integer, seed must be a positive integer\n');
    %error('ah');
    fn = NaN;
    FnVar = NaN;
    FnGrad = NaN;
    FnGradCov = NaN;
else
    
    %time after which orders of each product no longer arrive; i.e. Monday night
    %orders stop coming at 3 AM tuesday morning (t=27), etc.
    limits=[27;27;51;51;75;75;99;99;123;123;144;144;168];
    C=100;                                                      %Capacity
    b=x; %[100;100;100;100;100;100;100;100;100;100;100;100;100];%booking limits
    rates=(1/168)*[30,25,25,25,35,35,25,25,15,25,10,10,15];     %rates of arrival of orders
    price=100;                                                  %Fare paid/rate
    y=168;                                                      %Hours before t=0 at which reservations start to be taken.
    %i.e. if y=a, implies that orders start arriving at time -a.
    nProducts=13;
    
    %Incidence matrix as defined in problem statement: Rows are the resources
    %available, Columns are each "product". Entry ij=1 if product j uses
    %resource i, 0 otherwise.
    A=[ 1,1,0,0,0,0,0,0,0,0,0,0,0;
        0,1,1,1,0,0,0,0,0,0,0,0,0;
        0,0,0,1,1,1,0,0,0,0,0,0,0;
        0,0,0,0,0,1,1,1,0,0,0,0,0;
        0,0,0,0,0,0,0,1,1,1,0,0,0;
        0,0,0,0,0,0,0,0,0,1,1,1,0;
        0,0,0,0,0,0,0,0,0,0,0,1,1];
    
    
    TotalRevenue=0;
    %vector of next arrival times per product
    Arrival=zeros(nProducts,1)-y;
    
    arrBound = 10*round(168*sum(rates)); % upper bound on the number of arrivals over the time period
    aTime=zeros(nProducts,arrBound);
    a=ones(nProducts,1); % index of which arrival time to use next, for each product
    
    % Generate a new stream for random numbers
    ArrivalStream = RandStream.create('mrg32k3a');
    
    % Set the substream to the "seed"
    ArrivalStream.Substream = seed;
    
    % Generate random X
    OldStream = RandStream.setGlobalStream(ArrivalStream);
    
    for i = 1:nProducts
        aTime(i,:) = exprnd(1/rates(i), 1, arrBound);
    end
    
    RandStream.setGlobalStream(OldStream); %Return to old stream
    
    %Generate first arrivals.
    for i = 1:nProducts
        Arrival(i)=Arrival(i) + aTime(i, a(i));
        a(i)=a(i)+1;
    end
    
    minTime=0;
    while minTime <= runlength
        %Find next order received (minimum time) and make sure it is plausible
        %i.e. order time<limit for this product.
        minTime=runlength+1;
        for i=1:nProducts
            if(Arrival(i)<minTime && Arrival(i)<=limits(i))
                minTime=Arrival(i);
                minIndex=i;             %track product of next order.
            end
        end
        
        if minTime > runlength % Simulation time over, exit loop
            break
        end
        
        % Given next order, if b(i) > 0, accept it and update booking limits,
        % otherwise  reject the order.
        if b(minIndex)>0
            TotalRevenue=TotalRevenue+sum(price*A(:,minIndex));
            %Update booking limits: if two products share resources, decrease
            %booking lim
            for j=1:13
                if(A(:,minIndex)'*A(:,j)>=1)
                    if b(j)~=0
                        b(j)=b(j)-1;
                    end
                end
            end
        end
        %update time of next arrival for order received.
        Arrival(minIndex)=Arrival(minIndex)+aTime(minIndex, a(minIndex));
        a(minIndex)=a(minIndex)+1;
    end
    
    fn= TotalRevenue;
    FnVar=NaN;
    
end
end
