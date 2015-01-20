% main function, multi-frame blind deconvolution
% 
% call all the multi-frame blind deconvolution shell
% setup
localsetup;
% option.method = 'gaussian';
option.method = 'cg';
% generate psf
multiFilt = betterImRead; % 100 speckle samples
% sample psf
multiFilt_ds = multiFilt(1:1:50); 
numFrame = length(multiFilt_ds);
% generate multi frame with controlable noise
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option);
% mbd framework
start = multiFrame{1}; multiFrame = multiFrame(2:end);

% start = zeros(size(multiFrame{1})); % avarage of all frames
% for k = 1 : numFrame
%     start = start + multiFrame{k};
% end
% start = start./numFrame;

%
iter = 10;
I = mbd(multiFrame, F, start, iter, nature, multiFilt_ds,option);
