function clippedIm = clip(im,tolUp,tolDown)

if ~exist('tolUp', 'var'), tolUp=1e300; end
if ~exist('tolDown', 'var'), tolDown=1e-7; end

clippedIm = im;
clippedIm(im < tolDown) = tolDown;
clippedIm(im > tolUp) = tolUp;


end