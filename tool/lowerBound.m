function f = lowerBound(f)
% regularize kernel f
% manuelly setting small pixel value < 1/20 max(pixel value) = 0
%
maxPixelValue   =   max(vec(f));
factor          =   1/10;
f(f < factor * maxPixelValue) = 0;

end