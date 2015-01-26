function x = feasible(x)
% for image, make sure pixel value scaled in [0,1]
%
% xmed = median(vec(x));     % median of image
% x    = x./xmed * 0.5;      % scale wrt. median to 0.5 
x    = clip(x,1,0);        % clip

end