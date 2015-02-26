function f = lowerBound(f)
% regularize kernel f
% manuelly setting small pixel value < 1/20 max(pixel value) = 0
%
% maxPixelValue   =   max(vec(f));
% factor          =   1/10;
% f(f < factor * maxPixelValue) = 0;

%
maxPixelValue   =   max(vec(f));

w_row    = gausswin(size(f,1));
w_col    = gausswin(size(f,2))'; 

w        = bsxfun(@times, w_row, w_col);
w        = w./5;
w        = 0.2 - w;
w(w<0.1) = 0.05;
w(w>0.18) = 0.50;
fRef = maxPixelValue * w;
f(f < fRef) = 0;

% keyboard
% figure, imagesc(f)

end