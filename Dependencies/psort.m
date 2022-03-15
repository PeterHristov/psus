function [xS,parS,mnS,vrS] = psort(name,pars,mn,vr,x,n)
%%% Distributional sorting with modified first- and second-order stochastic dominance
%%% with Copeland counting.
%%% This version works with location-scale distributions only.
%%%
%%% [xS,parS,mnS,vrS] = psort(name,pars,mn,vr,x,n) takes the name of the distributions to sort - 
%%% 'Normal', 'Uniform', 'Logistic', 'Laplace', 'Student' (case insensitive), parameters
%%% of these distributions, mean and variance of the distributions, input data that generated
%%% the distributions, and the number of quadrature points for second order stochastic dominance.
%%% If used as part of psus, the parameters, means and variances come from the probabilistic code.
%%% The function outputs the input data, parameters, means and variances of the sorted distributions.

%% Preprocess for speed
if any(vr < 0), error('Variance must be non-negative.'); end

sd = sqrt(vr);
[distU,~,rInd] = unique([mn,sd,pars],'rows','stable'); %Find unique distros 
mnU = distU(:,1);
sdU = distU(:,2);
parsU = distU(:,3:end);

len = size(distU,1);
flag = zeros(len);

%% Work
for i = 1:len
    for j = i+1:len
        flag(j,i) = mfosd(mnU(j),mnU(i),sdU(j),sdU(i));
        if flag(j,i) == 0
            flag(j,i) = sosd(name,parsU(j,:),parsU(i,:),...
                mnU(j),mnU(i),sdU(j),sdU(i),n);
        end
    end
end

%% Build score matrix
flag = tril(flag) - tril(flag)';
fRank = sum(flag,2);
fRank = fRank(rInd); %Assign scores to repeated values

%% Sort
[~,sortOrd] = sort(fRank,'descend');

xS = x(sortOrd,:);
mnS = mn(sortOrd);
vrS = vr(sortOrd);
parS = pars(sortOrd,:);

end
