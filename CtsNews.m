%   ***************************************
%   *** Code written by German Gutierrez***
%   ***         gg92@cornell.edu        ***
%   ***                                 ***
%   *** Updated by Shane Henderson to   ***
%   *** use standard calling and random ***
%   *** number streams                  ***
%   ***************************************

% Last updated June 11, 2011


%Returns Mean of Profit and derivative. Arbitrarily return RH derivative

function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = CtsNews(x, runlength, seed, ~)
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = CtsNews(x, runlength, seed, other)
% x is the quantity to buy
% runlength is the number of days of demand to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used

constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;

if (x < 0) || (runlength <= 0) || (runlength ~= round(runlength)) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('x should be >= 0, runlength should be positive integer, seed must be a positive integer\n');
    fn = NaN;
    FnVar = NaN;
    FnGrad = NaN;
    FnGradCov = NaN;
else
    cost=5;
    sellPrice=9;
    salvage=1;
    alpha=2;
    beta=20;
    
    % Generate a new stream for random numbers
    OurStream = RandStream.create('mrg32k3a');
    
    % Set the substream to the "seed"
    OurStream.Substream = seed;
    
    % Generate demands
    OldStream = RandStream.setGlobalStream(OurStream);
    %OldStream = RandStream.setDefaultStream(OurStream); %versions 2010 and earlier
    demand = ((1-rand(runlength, 1)).^(-1/beta)-1).^(1/alpha);
    RandStream.setGlobalStream(OldStream); % Restore previous stream
    %RandStream.setDefaultStream(OldStream); %versions 2010 and earlier
    
    % Compute daily profit
    PerPeriodCost = cost * x;
    sales = min(demand, x);
    revenueSales = sellPrice * sales;
    revenueSalvage = salvage * (x - sales);
    profit = revenueSales + revenueSalvage - PerPeriodCost;
    RHDerivProfit = sellPrice * (demand > x) + salvage * (demand < x) - cost;
    fn = mean(profit);
    FnVar = var(profit) / runlength;
    FnGrad = mean(RHDerivProfit);
    FnGradCov = var(RHDerivProfit) / runlength;
end