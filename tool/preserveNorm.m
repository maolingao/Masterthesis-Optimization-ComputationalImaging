function f = preserveNorm(f)
% for PSF, 1. preserve energy norm;
%          2.make sure all element in f be positive.
%
f = clip(f, 1, 0);                      % clip btw [0,1]
f = f/(max(abs(vec(f))) + eps);         % scale down
f = f ./ (sum(vec(f)) + eps);           % preserve energy norm
%
end
%
%{
% unittest
% should give a energy-preserved f with sum of all elements in f = 1
f = rand(20);
f = preserveNorm(f);
sum(vec(f))
%}