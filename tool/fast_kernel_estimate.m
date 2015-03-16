function kernel_coarse = fast_kernel_estimate(X, gtImg, blurImg, option)
% X is the convolution operator class : conv2gradMat. It uses the image
% edge-taperred by [betterEdgeTaper.m] to generate itself.
% So the images to be convolved or correlated should also be edge-taperred 
% by [betterEdgeTaper.m]. --> suppress Gibbs effect.


% gtImg                           =   betterEdgeTaper(gtImg,option);                          % edge taper initial guess of g.t.
% blurImg                         =   betterEdgeTaper(blurImg,option);                        % edge taper initial guess of blurImg.

blurImgGrad                     =   cell(1,2);
[blurImgGrad{1},blurImgGrad{2}] =   gradient(blurImg);
gtImgGrad                       =   cell(1,2);
[gtImgGrad{1},gtImgGrad{2}]     =   gradient(gtImg);

kernel_coarse = (X' * blurImgGrad) ./ (X' * gtImgGrad + eps);
kernel_coarse = (X' * gtImgGrad + eps);

end