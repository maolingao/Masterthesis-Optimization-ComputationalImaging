% main function, multi-frame blind deconvolution
% 
% call all the multi-frame blind deconvolution shell
% setup
localsetup;
% option.method = 'gaussian';
% option.method = 'cg';
option.method = 'pncg';
% generate psf
multiFilt = betterImRead; % 100 speckle samples
% sample psf
multiFilt_ds = multiFilt(1:1:20); 
numFrame = length(multiFilt_ds);
% generate multi frame with controlable noise
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option);
% mbd framework
% --- use the first frame as start estimate for ground truth --- 
% start = multiFrame{1}; multiFrame = multiFrame(2:end);
%
% --- use the average of first 10 frames as start estimate for ground truth --- 
numFrameStart = 5;
start = zeros(size(multiFrame{1})); % avarage of all frames
for k = 1 : numFrameStart
    start = start + multiFrame{k};
end
start = start./numFrameStart; multiFrame = multiFrame((numFrameStart+1):end);
multiFilt_ds = multiFilt_ds((numFrameStart+1):end);

%
iter = 20;
I = mbd(multiFrame, F, start, iter, nature, multiFilt_ds,option);
