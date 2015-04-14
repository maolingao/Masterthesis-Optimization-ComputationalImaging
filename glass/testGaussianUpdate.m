% test gaussian
localsetup;

multiFilt   =   betterImRead; % 201 speckle samples
f           =   multiFilt{1};
f           =   f./sum(f(:));
sf          =   size(f);
%^^^^^^^
load('satel.mat');
% X           =   imread('cameraman.tif');
%^^^^^^^
x           =   im2double(X);  nature   =   x;
sx          =   size(x);
F           =   conv2MatOp(f,sx,'same');

im          =   F*x;
SNR         =   5;
start       =   ones(size(im));
im          =   addnoise( im, SNR, nature);
% figure, imagesc(im), colormap gray, axis image off

iter        =   50;
eta         =   0.1;

[gaussian_dI,errs,rerrs] = deconv_gaussian(F,im,iter,nature,start,eta,option);