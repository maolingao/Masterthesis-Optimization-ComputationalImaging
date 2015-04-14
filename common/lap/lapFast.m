function Lm = lapFast(m)
% fast calculate laplacian by shift-add 
% with boundary condition zero padding
% L = Lp - Ln
% L = [ 0 -1 0;
%       -1 4 -1;
%       0 -1 0];

Lm_neg = lapNegFast(m);
Lm_pos = 4*m;
Lm = Lm_pos - Lm_neg;
end