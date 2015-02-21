function im_tapered = betterEdgeTaper(im,option)
% edge tapering according to window type
%
if nargin < 2
    option.win = 'gaussian';
end
switch option.win
    case 'gaussian'
        w_row = gausswin(size(im,1));
        w_col = gausswin(size(im,2))'; 
      
    case 'barthann'
        w_row = barthannwin(size(im,1));
        w_col = barthannwin(size(im,2))'; 
        
    otherwise
        print '[betterEdgeTaper.m]: win can be either 'gaussian' or 'barthann''
end

w = bsxfun(@times, w_row, w_col);
% figure, imagesc(w), colormap gray, axis image
im_tapered = im .* w;

% figure, imagesc(im_tapered), colormap gray, axis image

end