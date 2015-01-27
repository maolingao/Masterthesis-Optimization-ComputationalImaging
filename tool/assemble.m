function [x_as] = assemble(x, x_1, alpha)
% assemble solution of 2 linear equations splitted in split.m
x_as = alpha .* x_1 + x;
end