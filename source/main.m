function convIm = main(shape)
%% generate blurred image
[nature,kernel] = initialize; % <-- guess ground truth
% [kernel,nature] = initialize; % <-- guess kernel

startup;
% option setup
option.version = 'FH';
option.figPath = '/is/ei/mgao/Documents/thesis/article/figure/lucy_regularization';

% call class
if nargin == 0
    shape = 'same';
end

xsize = size(nature);
F = conv2MatOp(kernel,xsize,shape);

% convolution
convIm = F*nature;
%% add noise
% SNR = 10;
% % convIm = addnoise(convIm, SNR, nature);
setupConv;
%% setups
iter =  5;
% initial guess
ci = 1; start = ci*(F'*convIm)+0*randn(xsize); start = start./sum(start(:)); % nfactor
% start = nature;
tol = 1e-20; 
eta = 0;
% ### call pncg ###
H = hessianMatrix(eye(size(F'*convIm)));
pncg_dI = deconv_pncg(F,convIm,nature,H,iter,start,tol,eta,option);
% ### call cg ###
cg_dI = deconv_cg(F,convIm,nature,iter,start,tol,eta,option);
% ### call lucy ###
lucy_dI = deconv_rl(F,convIm,iter,nature,start,eta,option);
% ### call gaussian ###
gaussian_dI = deconv_gaussian(F,convIm,iter,nature,start,eta,option);

% plot
saveResultFigure;
end
%%
function [im,f] = initialize
    % initiallize
    % ground truth
%     X = im2double(imread('cameraman.tif'));
    load('satel.mat');
    im = im2double(X);
    im = im./max(vec(im));
    % filter
%     f = fspecial('motion', 15, 30);
    load('AtmosphericBlur30.mat');
    f = PSF;
    scale = 31/length(f); f = imresize(f,scale);
    f = f/sum(f(:));
end


