function [flag,dX] = mfosd(mn1,mn2,sd1,sd2)
%%% Modified first order stochastic dominance for location-scale
%%% distributions.
%%% [flag,dX] = mfosd(mn1,mn2,sd1,sd2) - takes in the mean and standard deviations
%%% of two location-scale distributions to be compared and returns a flag = 1 if
%%% the first distribution is dominant, flag = -1 if the second distribution is 
%%% dominant, or flag = 0 if the two distributions are non-dominated to the first order.
%%% dX contains the differences that gave rise to the ranking.
 
%% Work
m = mn2-mn1;
s = 3*(sd2-sd1);

dX = [m-s,m];

if all(dX > 0), flag = -1;
elseif all(dX < 0), flag = 1;
else, flag = 0;
end