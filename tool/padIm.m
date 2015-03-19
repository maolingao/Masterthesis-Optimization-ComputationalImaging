function imPad = padIm(im,option)

if nargin < 2
    option.edgeTaper = 'blurExpansion';
end

xs = option.F.xsize;
fs = option.F.fsize;
ys = xs;

switch option.edgeTaper
    case 'taperEdge'   
        imsize = size(im);
        imPad = zeros(ys);
        imPad(1:imsize(1),1:imsize(2)) = im;
        imPad = circshift(imPad,ceil((ys - imsize)./2));  
        if sum(ys - imsize) ~= 0
            PSF = fspecial('gaussian',4*(ys - imsize),10);
            imPad = edgetaper(imPad,PSF);
        end    
        
    case 'blurExpansion'
%         imsize = size(im);
        f = fspecial('gaussian',min(xs,fs),10);
        F_blur = conv2MatOp(f,size(im),'same');
        tile = F_blur*im;
        
        imPad = [im fliplr(tile); flipud(tile) rot90(tile,2)];
%         figure, imagesc(imPad)
        imPad = circshift(imPad,ceil(size(tile)./2));  
%         figure, imagesc(imPad)      

    case 'taperEdgeBlurExpansion'
        
        imsize = size(im);
        f = fspecial('gaussian',4*min(xs,fs),10);
        imEdgeTaper = edgetaper(im,f);
        
        F_blur = conv2MatOp(f,imsize,'same');
        tile = F_blur*imEdgeTaper;
        
        imPad = [imEdgeTaper fliplr(tile); flipud(tile) rot90(tile,2)];
%         figure, imagesc(imPad)
        imPad = circshift(imPad,ceil(size(tile)./2));  
%         figure, imagesc(imPad)      

    case 'duplicateImage'
        imPad = [im fliplr(im); flipud(im) rot90(im,2)];
        imPad = circshift(imPad,ceil(size(im)./2));   

    otherwise
        error('check padIm')
end
        
% imPad = im;

end