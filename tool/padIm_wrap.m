function im_edgeTaper = padIm_wrap(im,option)
% edge taper images


im_edgeTaper = cellfun(@(x)padIm(x,option), im, 'UniformOutput', false);

end