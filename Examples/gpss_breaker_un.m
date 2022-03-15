function [y,u] = gpss_breaker_un(x,evals,mode)
% Functions specifically created and updated to try and break the GPSS
% algorithm to test its robustness. Based on Zhou (1998).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phi's
phi(:,1) = mvnpdf(x,0.04,[0.03,0.02;0.02,.025]);
phi(:,2) = mvnpdf(x,[0.98,0.7],[0.02,0;0,.003]);
phi(:,3) = mvnpdf(x,[0.75,0.85],[0.01,-0.015;-0.015,.03]);
phi(:,4) = mvnpdf(x,[0.71,0.32],[0.002,0;0,.002]);
phi(:,5) = mvnpdf(x,[0.33,0.83],[0.02,-0.01;-0.01,.01]);
phi(:,6) = mvnpdf(x,[0.43,0.73],[0.005,0;0,.005]);
phi(:,7) = mvnpdf(x,[0.23,0.93],[0.005,0;0,.005]);
phi(:,8) = mvnpdf(x,[1,0],[0.008,0;0,0.008]);
phi(:,9) = mvnpdf(x,[0.12,0.57],[0.005,0;0,.005]);

%% W
w(:,1) = 0.95;
w(:,2) = 0.278;
w(:,3) = 0.4167;
w(:,4) = 0.111;
w(:,5) = 0.85*0.55;
w(:,6) = 0.85*0.08;
w(:,7) = 0.85*0.09;
w(:,8) = 0.302;
w(:,9) = 0.236;

sumW = sum(w);
w = w/sumW;
%% Sum
m = sum(phi*diag(w),2);

%% Uncertainty
if strcmp(mode,'simple')
    u = evals^-1*ones(size(x,1),1); %Variance of a normal deviate
    y = m - sqrt(u); %To enable early subdomain discovery
else
    cv = cv_uncert(x);
    u = ( (50/evals) *cv .* m).^2;
    y = m - sqrt(u); %Almost exact convergence
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