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

multiFrame = cell(1,numFrame);
[nature,~] = initialize;
xsize = size(nature);
shape = 'same';
% figure(33), clf

for i = 1:numFrame + 1

    if i == numFrame+1
        break;
    end
    
%     multiFilt{i} = kernel;    % ---> change for rotateKernel
    multiFilt{i} = im2double(multiFilt{i});    
    scale = 21/length(multiFilt{i}); multiFilt{i} = imresize(multiFilt{i},scale);
    multiFilt{i} = multiFilt{i}./sum(vec(multiFilt{i}));
    kernel = multiFilt{i};
    % save kernel image  
    %{  
    fkernel = figure; set(fkernel,'visible','off');
    imagesc(clip(kernel,1,0)); axis image, colormap(gray)
    filename = sprintf('kernel_%d_frame',i);
    filename = fullfile(option.figPath,filename);
    print(gcf, '-depsc2', filename)
    %}
    % ###############
    % check
%     kernel = imresize(kernel,20/length(kernel));
    % ###############
    kernel = kernel./sum(vec(kernel));
    multiFilt{i} = kernel;
    
    F = conv2MatOp(kernel,xsize,shape); % <-- convolution matrix
    convIm = F*nature; % <-- convolvution operation
    % add noise
%     SNR = 5;
%     convIm = addnoise(convIm,SNR,nature);
    figure(333), imagesc(convIm), colormap gray, axis off image, 
    if exist('SNR','var')
        title(sprintf('SNR = %ddB',SNR));
    else
        title(sprintf('blurry image'));
    end
    % ################
    % check
    %{
    eps = 1e-7;
    [convMtx,convMtxTranpose] = buildF(shape,kernel,nature);
% %     x_ch = reshape(((convMtx'*convMtx + eps)\(convMtx'*vec(convIm) + eps)), xsize);   % <--- F
% %     diff = nature - x_ch;                                                             % <--- F
    fsize = size(kernel);                                                                % <--- X
    kernel_ch = reshape(((convMtx'*convMtx + eps)\(convMtx'*vec(convIm) )),fsize);  % <--- X
    diff = kernel - kernel_ch;                                                           % <--- X
    fprintf('kernel difference:%d\n',max(vec(abs(diff))))
    if i == 2
        save convMtx.mat convMtx
    end
    
% %     convMtx_ch = (vec(convIm)*transpose(vec(nature)))*(inv(vec(nature)*transpose(vec(nature)))); % <--- F
    convMtx_ch =  (kernel(:)*kernel(:)'+ eps) \ (convMtx'*vec(convIm)*kernel(:)'); % <--- X
    diff = convMtx'*convMtx - convMtx_ch;                             % <--- X
    fprintf('convMtx^T*convMtx difference:%d\n',max(vec(abs(diff))))
    %}
    % ################
    multiFrame{i} = convIm;
    % ------------ show kernel and blurry image ------------ %
    %{
    figure(33), subplot(numFrame+1,2,i*2-1),imshow(imresize(kernel/max(kernel(:)),20))
    title(sprintf('kernel %d/%d',i,numFrame))
    subplot(numFrame+1,2,i*2),imshow(convIm/max(convIm(:)))
    title(sprintf('image %d/%d',i,numFrame))
    drawnow
    %}
    
%     kernel = rotatekernel(kernel,theta); % shift kernel every 5 degree  % ---> change for rotateKernel
    
    
%     [~,kernel] = initialize(len,i+angle); % motion kernel in different directions
%     kernel = circshift(kernel,[0,1]);
%     figure(14)
%     imshow(imresize(kernel/max(kernel(:)),20))
%     title(sprintf('rotated kernel'))
%     drawnow
    

end

% ####### figure ######
setupConv;
end

function [im,f] = initialize
    % initiallize
%     ground truth
    X = im2double(imread('cameraman.tif'));
%     load('satel.mat');
    im = im2double(X);
    im = im./max(vec(im));
    im = im./10^(3);  % <---- scale down image
    % filter
%     f = fspecial('motion', 15, 30);
    load('AtmosphericBlur30.mat');
    f = PSF;
    scale = 31/length(f); f = imresize(f,scale);
    f = f/sum(f(:));
end