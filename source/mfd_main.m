% main function
% 
% call all the multi-frame deconvolution shell
% setup
localsetup;
figure(200), clf
% option.method = [{'cg'},{'pncg'},{'gaussian'},{'rl'}];
% option.method = [{'cg'},{'pncg'},{'gaussian'}];
% option.method = [{'cg'},{'pncg'}];
% option.method = [{'pncg'},{'gaussian'},{'rl'}];
% option.method = [{'pncg'},{'gaussian'}];
% option.method = [{'gaussian'}];
option.method = {'pncg'};
% option.method = [{'cg'}];
option.version = 'FH';      % which model of pncg to use
option.target = 'kernel';   % what object to guess
option.color = 'dre';       % set color and linestyle of hData
option.LineStyle = '-';     % set color and linestyle of hData
% generate psf
multiFilt = betterImRead; % speckle blur
% sample psf
multiFilt_ds = multiFilt(1:1:20); 
numFrame = length(multiFilt_ds);
% generate multi frame with controlable noise
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option); % ---> change for rotateKernel
% mfd framework
% start = multiFrame{1}; multiFrame = multiFrame(2:end);
%
iter = 10; %length(nature); %numel(multiFilt_ds{1}); % 
I = mfd(multiFrame,multiFilt_ds,F,iter,nature,option);
% close all hidden figures
% close all hidden