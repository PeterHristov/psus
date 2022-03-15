function [flag,dF,intF1,intF2,x] = sosd(name,pars1,pars2,mn1,mn2,sd1,sd2,n)
%%% Second order stochastic dominance.
%%% [flag,dF,intF1,intF2,x] = sosd(name,pars1,pars2,mn1,mn2,sd1,sd2,n) compares
%%% two distributions via second order stochastic dominance, by taking
%%% in the name of the two-parameter distributions and their means and variances.
%%% The number of quadrature points is specified through the input n.
%%% The function returns the dominance flag - 1, if the first distribution is dominant
%%% -1, if the second distribution is dominant, and 0 if the two distributions are
%%% non-dominated. To inspect the result, also returns dF, whcih contains the differences 
%%% between integrals intF1 and intF2, and the quadrature points in the vector x.
%%%
%%% The function currently supports only location-scale distributions in line with the use
%%% of psus.

xmin = min(mn1 - 5*sd1,mn2 - 5*sd2);
xmax = max(mn1 + 5*sd1,mn2 + 5*sd2);
x = linspace(xmin,xmax,n);

intF1 = zeros(n,1);
intF2 = zeros(n,1);

[cdf1,cdf2] = getCDF(name,pars1,pars2); %Generate appropriate CDFs

for i = 1:n-1
    intF1(i+1,1) = trapz(x(1:i+1),cdf1(x(1:i+1)));
    intF2(i+1,1) = trapz(x(1:i+1),cdf2(x(1:i+1)));
end

dF = intF1 - intF2;
sdF = sum(dF);

if sdF > 0, flag = -1; %Second one dominant
elseif sdF < 0, flag = 1; %First one dominant
else, flag = 0; %Indistinguishable => identical
end

end

function [cdf1,cdf2] = getCDF(name,pars1,pars2)
outDists = {'Normal'; 'Uniform'; 'Logistic'; 'Laplace'; 'Student'};

name = outDists{strncmpi(name,outDists,4)};

switch name
    case 'Normal'
        cdf1 = @(x) normcdf(x,pars1(1),pars1(2));
        cdf2 = @(x) normcdf(x,pars2(1),pars2(2));
    case 'Uniform'
        cdf1 = @(x) unifcdf(x,pars1(1),pars1(2));
        cdf2 = @(x) unifcdf(x,pars2(1),pars2(2));
    case 'Logistic'
        cdf1 = @(x) cdf('logistic',x,pars1(1),pars1(2));
        cdf2 = @(x) cdf('logistic',x,pars2(1),pars2(2));
    case 'Laplace'
        cdf1 = @(x) lapcdf(x,pars1(1),pars1(2));
        cdf2 = @(x) lapcdf(x,pars2(1),pars2(2));
    case 'Student'
        cdf1 = @(x) tcdf(x-pars1(1),pars1(2));
        cdf2 = @(x) tcdf(x-pars2(1),pars2(2));
    otherwise
        error('Unrecognized output distribution.')
end

end

function p = lapcdf(lev,m,s)

if lev <= m
    p = 0.5*exp( (lev-m)./s );
else
    p = 1 - 0.5*exp( (m-lev)./s );
end

end
