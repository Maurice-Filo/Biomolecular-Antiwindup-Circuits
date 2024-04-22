function h = HillFunction(x, alpha, kappa, n)
h = alpha * x.^n ./ (x.^n + kappa^n);
end

