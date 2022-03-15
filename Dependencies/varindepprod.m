function varprod = varindepprod(mu,var)
%%% Variance of the product of two independent distributions.
%%% Inputs:
%%% 	mu -  a vector of expectations
%%% 	var - a vector of variances

n = length(mu); %Number of RV

if n == 1, varprod = var;
elseif n == 2
    varprod = prod(var) + var(1)*mu(2)^2 + var(2)*mu(1)^2;
else
    v = varindepprod(mu(1:n-1),var(1:n-1));
    varprod = var(n)*v + v*mu(n)^2 + var(n)*prod(mu(1:n-1).^2);
end

end