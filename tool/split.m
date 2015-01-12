function [b_orth, alpha] = split(b,basis)
% split A*x = b into 2 linear problems: 1. A*x = alpha*basis
%                                       2. A*x = b_orth
% namely a linear equation in the basis direction and its orthogonal one.

alpha = proj(b,basis);
b_orth = b - alpha * basis;

end