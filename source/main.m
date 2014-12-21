function convIm = main(shape)
%% generate blurred image
[nature,kernel] = initialize; % <-- guess ground truth
% [kernel,nature] = initialize; % <-- guess kernel

startup;

% call class
if nargin == 0
    shape = 'full';
end

xsize = size(nature);
F = conv2MatOp(kernel,xsize,shape);

% convolution
convIm = F*nature;
%% add noise
%% setups
iter =  100;
% initial guess
ci = 1; start = ci*(F'*convIm)+0*randn(xsize); start = start./sum(start(:)); % nfactor
tol = 1e-20; 
eta = 0.5;
% ### call pncg ###
H = hessianMatrix(eye(size(F'*convIm)));
pncg_dI = deconv_pncg(F,convIm,nature,H,iter,start,tol,eta);
% ### call cg ###
cg_dI = deconv_cg(F,convIm,nature,iter,start,tol,eta);
% ### call lucy ###
lucy_dI = deconv_rl(F,convIm,iter,nature,start,eta);
% ### call gaussian ###
gaussian_dI = deconv_gaussian(F,convIm,iter,nature,start,eta);

end
%%
function [im,f] = initialize
    % initiallize
    % ground truth
    X = im2double(imread('cameraman.tif'));
    im = im2double(X);
    im = im./max(vec(im));
    % filter
    f = fspecial('motion', 15, 30);
    scale = 31/length(f); f = imresize(f,scale);
    f = f/sum(f(:));
end


