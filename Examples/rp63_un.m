function [y,u] = rp63_un(x,evals)
%%% An d-dimensional function, presented at the 2019 Reliability Black-Box
%%% challenge. The function is vectorised - x is an n-by-d array and evals
%%% is the number of itereations from the solver.

m = -(0.1 * sum(x(:,2:end).^2,2) - x(:,1));

u = evals^-1*ones(size(x,1),1);

y = m - sqrt(u); 
end