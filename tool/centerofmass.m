function v = centerofmass(f, support);
% functionality: 
% support = 1 : pure    geometrical center of all mass points
%         = 0 : weigted geometrical center of all mass points
if ~exist('support', 'var') || isempty(support), support = 0; end
if support == 1, f(f>0) = 1; end
sf = size(f);
f  = f / sum(f(:));
i  = sum((1:sf(1))*f);
j  = sum((1:sf(2))*f');
v  = [i,j];
return
