function [multiFrame,multiFilt,F,nature] = generateMultiFrame(numFrame,multiFilt,option)
% build multi blurred image frame
% return blurry images and corresponding kernels

if nargin < 1
    numFrame = 3;
end
if nargin < 2
    multiFilt = cell(1,numFrame); % register for filters
end

if nargin < 3
    option.figPath = '/is/ei/mgao/figure2drag';
end
if isfield(option,'SNR')
    SNR     =   option.SNR;
end

multiFrame  =   cell(1,numFrame);
[nature,~]  =   initialize;
xsize       =   size(nature);
shape       =   'same';

for i = 1:numFrame + 1

    if i == numFrame+1
        break;
    end
    
    % ########## generate blurry frames ##########
    multiFilt{i}    =  im2double(   multiFilt{i} );    
    scale           =  21 / length( multiFilt{i} );                  % sample kernel, make a smaller size
    multiFilt{i}    =  imresize(    multiFilt{i}, scale );
    multiFilt{i}    =  multiFilt{i} ./ sum( vec( multiFilt{i} ));    % upper bound pixel value = 1
    kernel          =  multiFilt{i};
    
    F               =  conv2MatOp(kernel,xsize,shape);               % <-- convolution matrix
    convIm          =  F * nature;                                   % <-- convolvution operation
    % ########## add noise ##########
    if isfield(option,'SNR')
        convIm      =  addnoise( convIm, SNR, nature);
    end
    multiFrame{i}   =  convIm;
    % ########## show kernel and blurry image ##########
    %{
    figure(33), subplot(numFrame+1,2,i*2-1),imshow(imresize(kernel/max(kernel(:)),20))
    title(sprintf('kernel %d/%d',i,numFrame))
    subplot(numFrame+1,2,i*2),imshow(convIm/max(convIm(:)))
    title(sprintf('image %d/%d',i,numFrame))
    drawnow
    %}    

end

    figure(333), imagesc(convIm), colormap gray, axis off image;
    if exist('SNR','var')
        title(sprintf('SNR = %ddB',SNR));
    else
        title(sprintf('blurry image'));
    end

% ####### save setup figures ######
setupConv;
end

function [im,f] = initialize
    % initiallize
%     ground truth
%     X   =  im2double(imread('cameraman.tif'));
    load('satel.mat');
    im  =   im2double(X);
    im  =   im./max(vec(im));
    im  =   im./10^(3);                         % <---- scale down image
    % filter
%     f   =   fspecial('motion', 15, 30);
    load('AtmosphericBlur30.mat');
    f     =   PSF;
    scale =   31/length(f); 
    f     =   imresize(f,scale);
    f     =   f/sum(f(:));
end