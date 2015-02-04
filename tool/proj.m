function alpha = proj(b,basis)
% project b onto basis, return alpha which satisfies b_proj = alpha * basis

epsl = 1e-8;
alpha = (basis' * b) / (basis' * basis + epsl);

end