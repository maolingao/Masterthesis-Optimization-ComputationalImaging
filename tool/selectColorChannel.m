function videoFrameMono = selectColorChannel(videoFrame, channel)
% extract predefined color channel
% also valid for gray scale images

if ~exist('channel','var')
    channel = 1; % red 
end

videoFrameMono = cellfun(@(x)x(:,:,channel), videoFrame, 'UniformOutput', false);

end