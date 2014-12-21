function clippedIm = clip(im,tolUp,tolDown)

clippedIm = im;
clippedIm(im < tolDown) = tolDown;
clippedIm(im > tolUp) = tolUp;


end