function A = RandomRotation(D)
% produces a rotation matrix A in SO(D) (the special orthogonal
% group SO(D), or orthogonal matrices with unit determinant, drawn
% uniformly from the Haar measure. 
%
% The algorithm used is the subgroup algorithm as originally proposed by 
% 
% P. Diaconis & M. Shahshahani, "The subgroup algorithm for generating
% uniform random variables". Probability in the Engineering and
% Informational Sciences 1: 15?32 (1987)
%
% This is not a particularly efficient implementation.
%
% code by Philipp Hennig, November 2013

assert(D >= 2);

% induction start: uniform draw from D=2 Haar measure

t = 2*pi*rand();
A = [cos(t), sin(t); -sin(t), cos(t)];

% (to draw from the Haar measure on O(D), instead of SO(D), draw
% if rand() > 0.5
%     A = [cos(t), sin(t); -sin(t), cos(t)];
% else
%     A = [cos(t), sin(t); sin(t), -cos(t)];
% end

for d = 2:D-1 % induction step
    v = randn(d+1,1); 
    v = v ./ sqrt(v'*v); % draw on S_d the unit sphere
    e = [1; zeros(d,1)];
    x = (e - v) ./ sqrt((e-v)' * (e-v)); % random coset location of SO(d-1) in SO(d)
    
    D = [1, zeros(1,d); zeros(d,1), A]; 
    A = D - 2 * x * (x' * D);
end

A = -A; % fix determinant.