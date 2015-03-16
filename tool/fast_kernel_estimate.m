function kernel_coarse = fast_kernel_estimate(X, gtImg, blurImg, option)
% X is the convolution operator class : conv2gradMat. It uses the image
% edge-taperred by [betterEdgeTaper.m] to generate itself.
% So the images to be convolved or correlated should also be edge-taperred 
% by [betterEdgeTaper.m]. --> suppress Gibbs effect.
% input:
% X         : class conv2gradMat, for retrieving convolution information 
% gtImg     : edge-taperred ground truth image / ground truth estimation
% blurImg   : edge-taperred blurry image 
%
% output:
% kernel_coarse : estimated coarse kernel


if isprop(X,'GX')
    Pf = X.GX{1}.Pf;
    Px = X.GX{1}.Px;
    Py = X.GX{1}.Py;
else
    Pf = X.Pf;
    Px = X.Px;
    Py = X.Py;
end
    
% gradient images
blurImgGrad                     =   cell(1,2);
[blurImgGrad{1},blurImgGrad{2}] =   gradient(blurImg);
gtImgGrad                       =   cell(1,2);
[gtImgGrad{1},gtImgGrad{2}]     =   gradient(gtImg);

% pad images
blurImgGrad         =   cellfun(@(x)Py*x, blurImgGrad, 'UniformOutput', false);
gtImgGrad           =   cellfun(@(x)Px*x, gtImgGrad, 'UniformOutput', false);

% fft images
blurImgGradFT       =   cellfun(@fft2, blurImgGrad, 'UniformOutput', false);
gtImgGradFT         =   cellfun(@fft2, gtImgGrad, 'UniformOutput', false);
gtImgGradFTconj     =   cellfun(@conj, gtImgGradFT, 'UniformOutput', false);

% conj(fft2(x)) .* fft2(y)
numeratorCellarray  =   cellfun(@(x1,x2)x1.*x2, gtImgGradFTconj , blurImgGradFT, 'UniformOutput', false);
denominatorCellarray =  cellfun(@(x1,x2)x1.*x2, gtImgGradFTconj , gtImgGradFT, 'UniformOutput', false);

% sum up of 2 gradient images with different direction
numerator           =   numeratorCellarray{1} + numeratorCellarray{2};
denominator         =   denominatorCellarray{1} + denominatorCellarray{2};

% analytic kernel solution
kernel_coarse       =   Pf'*(real(ifft2(numerator ./ (denominator + eps))));

% positivity and energy preserving
kernel_coarse       =   kernel_coarse - min(vec(kernel_coarse));
kernel_coarse       =   preserveNorm(kernel_coarse);

end

% ^^^^^^^^^ unittest
%{
% should return a number close to zero
option.win          =   'barthann';
load('AtmosphericBlur30.mat');  f = PSF;
natureI = im2double(imread('cameraman.tif'));
natureI = natureI./max(vec(natureI));
X = conv2MatOp(natureI,size(f),'same');
frame   = X * f;
natureI4convmat     =   betterEdgeTaper(natureI,option);
frameEdgeTaperred   =   betterEdgeTaper(frame,  option);
kernel_coarse       =   fast_kernel_estimate(X, natureI4convmat, frameEdgeTaperred, option);
max(vec(kernel_coarse - f))
%}
% ^^^^^^^^^ 