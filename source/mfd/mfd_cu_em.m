function I = mfd_cu_em(multiFrame, F, start, iter, natureI, multiKernel, option)
% [pncg_dI, cg_dI, lucy_dI, gaussian_dI] = mbd(multiFrame, F, start, iterK, natureI, multiKernel)
% offline multi-frame blind deconvolution
if nargin < 7
    option.figPath = '/is/ei/mgao/figure2drag';
    option.method = 'gaussian';
    option.win = 'barthann';
    option.MEMLIM = 50;
    option.img = 'gradient';
    option.blind = 'nb';
end

startup;
figPath = option.figPath;

figure(200), clf

if F.xsize > F.fsize
    imagesize = F.xsize;
    fsize = F.fsize;
else
    imagesize = F.fsize;
    fsize = F.xsize;
end
shape = F.shape;
% ##### general setup #####
numFrame = numel(multiFrame);
for i = 1 : numFrame                            % register observation images
    fixed = natureI;                            % r.t. ground truth
    moving = multiFrame{i};
    subpixel = 1;
    [multiFrame{i}, output] = efficient_imregister(fixed, moving, subpixel);
    %
    multiKernel{i} =  center(multiKernel{i});   % center all kernel, only for error analysis
end
eta = option.eta;                                        % ### <--- regularization parameter
% kernel estimating
tol    =  option.tolK;
time  = 0;

natureI4convmat   =   betterEdgeTaper(natureI,option);
%%%%%%%%%%%%%

% iteration
methodsAmount = length(option.method);
switch option.mode
    case 'compPerFrame'
        loopCompPerFrame = numFrame;
        loopCompAllFrame = 1;
    case 'compAllFrame'
        loopCompPerFrame = 1;
        loopCompAllFrame = numFrame;
    otherwise
        display('in [mfd_cu.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
end
        
for j = 1 : loopCompPerFrame
for i = 1 : methodsAmount
    %%
        f_pncg = figure(114); clf(f_pncg); set(f_pncg,'visible','on')
        for jj = 1 : loopCompAllFrame            
            switch option.mode
                case 'compPerFrame'
                    k = j;
                case 'compAllFrame'
                    k = jj;
            end
            % ##### special setup #####
            frame       =   clip(multiFrame{k},inf,0);
            natureK     =   multiKernel{k};                 % for error calculation
            
            % !!!!!!!! edge taper !!!!!!!!
            frameEdgeTaperred  =   betterEdgeTaper(frame,  option);
            switch option.img
                case 'normal'
                    % !!!!!!!! normal image !!!!!!!!
                    X                 =   conv2MatOp(im2double(natureI4convmat),fsize,shape);              % convMtx of nature
                    frame4estiKernel  =   frameEdgeTaperred;
                    startK            =   X'* frame ;                                                      % initial guess of kernel, b
                    % !!!!!!!! gradient image !!!!!!!!
                case 'gradient'
                    X                           =   conv2gradMat(im2double(natureI4convmat),fsize,shape);  % convMtx of gradient image of nature
                    frameGrad                   =   cell(1,2);
                    [frameGrad{1},frameGrad{2}] =   gradient(frameEdgeTaperred);
                    frame4estiKernel            =   frameGrad;
                    startK                      =   X'* frameGrad;                                         % initial guess of kernel, b
            end
            
            option.plotFlag     =   1;
            switch option.solver
                case 'gaussian'
                    [dI_kernel,data_k] = deconv_gaussian(X,frame4estiKernel,iter,natureK,startK,eta,option); % gaussian
                case 'lucy'
                    [dI_kernel,data_k] = deconv_rl(X,frame4estiKernel,iter,natureK,startK,eta,option); % richardson-lucy
            end
            dI_kernel         =   preserveNorm(dI_kernel);            % preserve energy norm of PSF
            % ###### psf residual curve ######
            data.errs(:,k) = data_k.errs;
            data.rerrs(:,k) = data_k.rerrs;
            data.time(:,k) = data_k.time;
            figure(200), set(gcf,'visible','on')
            subplot(121), loglog(1:length(data_k.errs),data_k.errs); hold on
            subplot(122), loglog(data_k.time,data_k.errs); hold on
            
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,dI_kernel,k,'pncg','cmp',figPath);
       
            % statitics 
            % all frame errors - kernel
            time               = timeLabel(time, data_k.time);
        end
            I = dI_kernel;
            
            
end 
switch option.solver
    case 'gaussian'
        save('data_nb_gaussian.mat','data');
    case 'lucy'
        save('data_nb_rl.mat','data');
end
%##### set, label and save figure #####
saveResultFigure;

end
end
