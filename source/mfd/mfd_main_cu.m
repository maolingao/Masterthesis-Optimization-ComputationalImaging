% main function, multi-frame blind deconvolution
% 
% call all the multi-frame blind deconvolution shell
% --------------------    setup   -------------------- %
localsetup;
% -------------------- readin PSF -------------------- %
% generate psf
multiFilt    = betterImRead; % 100 speckle samples
% load('/is/ei/mgao/Documents/thesis/Astro/simulation/motionblur/multiMotionBlur.mat');
% multiFilt = multiMotionBlur;
% sample psf
numFrame     = option.numFrame;
% multiFilt_ds = multiFilt(randperm(length(multiFilt),numFrame)); 
multiFilt_ds = multiFilt(1:numFrame); 
% multiFilt_ds = multiFilt(1:numFrame); 
% -------------- generate blurry frames -------------- %
% generate multi frame with controlable noise
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option);
% --------------------   mfd     -------------------- %
% mfd framework
start        =   nature; 
%
%
iter =  option.iter;
I    =  mfd_cu(multiFrame, F, start, iter, nature, multiFilt_ds, option);
option.solver = 'lucy';
I    =  mfd_cu_em(multiFrame, F, start, iter, nature, multiFilt_ds, option);
