% main function, multi-frame blind deconvolution
% 
% call all the multi-frame blind deconvolution shell
% --------------------    setup   -------------------- %
localsetup;
% -------------------- readin PSF -------------------- %
% generate psf
multiFilt    = betterImRead; % 100 speckle samples
% sample psf
numFrame     = option.numFrame;
multiFilt_ds = multiFilt(1:1:numFrame); 
% -------------- generate blurry frames -------------- %
% generate multi frame with controlable noise
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option);
% --------------------   mbd     -------------------- %
% mbd framework
% initial estimate for ground truth 
numFrame4Start = option.numFrame4Start;   % num of blurry frame for initial estimate
start        =   zeros(size(multiFrame{1}));     % avarage of those frames
for k = 1 : numFrame4Start
    start    =   start + multiFrame{k};
end
start        =   start./numFrame4Start; 
multiFrame   =   multiFrame((numFrame4Start+1):end);
multiFilt_ds =   multiFilt_ds((numFrame4Start+1):end);
%
%
iter =  option.iter;
I    =  mbd(multiFrame, F, start, iter, nature, multiFilt_ds, option);
