% main function, multi-frame blind deconvolution
% 
% call all the multi-frame blind deconvolution shell
% --------------------    setup   -------------------- %
localsetup;

% -------------------- readin video frames -------------------- %
% generate psf
path = '/is/ei/mgao/Documents/thesis/Astro/real_data/star';
videoFrame    = betterImRead(path); 
% frame amount
if option.numFrame ~= inf
    videoFrame = videoFrame(1:option.numFrame);
else
    option.numFrame = length(videoFrame);
end
% extract color channel
videoFrameMono = selectColorChannel(videoFrame, 1);
% scale images
[videoFrameMono] = scaling(videoFrameMono, 1e-3);

% --------------------   mbd     -------------------- %
% mbd framework
% initial estimate for ground truth 
numFrame4Start = option.numFrame4Start;   % num of blurry frame for initial estimate
start        =   zeros(size(videoFrameMono{1}));     % avarage of those frames
for k = 1 : numFrame4Start
    start    =   start + videoFrameMono{k};
end
start        =   start ./ (numFrame4Start + eps); 
videoFrameMono   =   videoFrameMono((numFrame4Start+1):end);
%
%
iter =  option.iter;
nature = videoFrameMono{1};  % for register all frames
PSFs = cell(option.numFrame - numFrame4Start,1);
PSFs = cellfun(@(x) eye(option.F.fsize), PSFs, 'UniformOutput', false);
I    =  mbd(videoFrameMono, option.F, start, iter, nature, PSFs, option);
