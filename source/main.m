function convIm = main(shape)
localsetup;

%% generate blurred image
[nature,kernel] = initialize; % <-- guess ground truth
% [kernel,nature] = initialize; % <-- guess kernel

startup;
% option setup
option.version = 'FH';
% option.figPath = '/is/ei/mgao/Documents/thesis/article/figure/lucy_regularization';
% option.color = 'ora';
option.LineStyle = '-';

% call class
if nargin == 0
    shape = 'same';
end

%% --- regularization figure setting ---
colors = {'dre','ora','blu','gra','mpg','tumblue','lightdre', 'lightora', 'lightblu','lightmpg'};
LineStyles = {'-','--','-.','--','-','--','-.','--','-',':'};
figure(34), clf
for i=1:10
option.color = colors{i};
option.LineStyle = LineStyles{i};

%% 
xsize = size(nature);
F = conv2MatOp(kernel,xsize,shape);
% convolution
convIm = F*nature;
%% add noise
SNR = 10;
[convIm, option.noiseVar] = addnoise(convIm, SNR, nature);
setupConv;
%% setups
iter =  50;
% initial guess
ci = 1; start = ci*(F'*convIm)+0*randn(xsize); start = start./max(start(:)); % nfactor, if guessing Kernel!
% start = nature;
tol = 1e-20; 
eta = (i - 1)./(i - 1 + eps) * 10^(-7 + i);
% eta = 0;
% ### call pncg ###
% H = hessianMatrix(eye(size(F'*convIm)));
% pncg_dI = deconv_pncg(F,convIm,nature,H,iter,start,tol,eta,option);
% ### call cg ###
% cg_dI = deconv_cg(F,convIm,nature,iter,start,tol,eta,option);
% ### call lucy ###
% keyboard
lucy_dI = deconv_rl(F,convIm,iter,nature,start,eta,option);
% ### call gaussian ###
% gaussian_dI = deconv_gaussian(F,convIm,iter,nature,start,eta,option);

end
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


