function [f,u] = branin_uncert(x,evals,mode)
%x is a 1x2 vector 
%
if nargin < 2, evals = 1; end

%% Input scaling
u1 = x(:,1); %Uniform scales
u2 = x(:,2);

x1 = 15*u1-5;
x2 = 15*u2;

%% Mean
a = 1;
b = 5.1/4/pi^2;
c = 5/pi;
r = 6;
s = 10;
t = 1/(8*pi);

m = (a*(x2-b*x1.^2+c*x1-r).^2+s*(1-t)*cos(x1)+s)+5*x1; %True function

%% Uncertainty
if strcmpi(mode,'const') %Variance depends only on iterations
    %evals = min(evals,100); %Prevent negative variance
    %u = 50./evals-.5 * ones(size(x,1),1); %For example
    cv = cv_uncert(x);
else
    u = sqrt(exp(-0.01*evals) * abs(m)); %SD of a normal deviate
end
%% Mean of observation
%scales = log(abs(m)+1) / 2; %Original c
scales = min(log(abs(m)+1) / 2, 1.8); %Modified c for UC

if strcmpi(mode,'pos')
    f = m + scales.*u;
elseif strcmpi(mode,'neg')
    f = m - scales.*u;
else
    u = (cv .* m);
    f = m; %Almost exact convergence
end

end

function cv = cv_uncert(x)
cv = zeros(size(x,1),1);

for i = 1:size(x,1)
    nX = norm(x(i,:));
    
    if nX<0.2
        cv(i,1) = 0.0028;
    elseif nX<0.4 && nX>0.2
        cv(i,1) = 0.0049;
    elseif nX<0.6 && nX>0.4
        cv(i,1) = 0.0946;
    elseif nX<0.8 && nX>0.6
        cv(i,1) = 0.0151;
    elseif nX<1 && nX>0.8
        cv(i,1) = 0.0023;
    elseif nX<1.2 && nX>1
        cv(i,1) = 0.0102;
    elseif nX<1.4 && nX>1.2
        cv(i,1) = 0.0291;
    elseif nX<1.6 && nX>1.4
        cv(i,1) = 0.0420;
    elseif nX<1.8 && nX>1.6
        cv(i,1) = 0.0464;
    elseif nX<2 && nX>1.8
        cv(i,1) = 0.0531;
    end
end
end