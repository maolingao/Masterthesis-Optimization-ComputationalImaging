function [b,c] = swap(b,c)
% swap the values of 2 variables
c_reg = c;
c = b;
b = c_reg;
end