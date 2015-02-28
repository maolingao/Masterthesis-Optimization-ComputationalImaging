function I = mbd(multiFrame, F, start, iterK, natureI, multiKernel, option)
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

f10 = figure(10); clf(10), set(f10,'visible','off')
f11 = figure(11); clf(11), set(f11,'visible','off')
f12 = figure(12); clf(12), set(f12,'visible','off')
f13 = figure(13); clf(13), set(f13,'visible','off')
fclk = figure(14); clf(fclk), set(fclk,'visible','on')
fstp = figure(15); clf(fstp), set(fstp,'visible','on')

if F.xsize > F.fsize
    imagesize = F.xsize;
    fsize = F.fsize;
else
    imagesize = F.fsize;
    fsize = F.xsize;
end
shape = F.shape;

errs_allframes_pncg = []; rerrs_allframes_pncg = [];
errs_allframes_pncgK = []; rerrs_allframes_pncgK = [];
errs_allframes_cg = []; rerrs_allframes_cg = [];
errs_allframes_gaussian = []; rerrs_allframes_gaussian = [];
errs_allframes_rl = []; rerrs_allframes_rl = [];
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
tolK   =  -inf;
scaler =  1e0;
HK     =  hessianMatrix(eye(fsize)*scaler);
start  =  clip(start,inf,0);
start(start<1e-7) = 0;
cutLine = nan;
% start  =  start./max(vec(start)) * 1e-3;
%%%%%%%%%%%%%
fixed    =  natureI;                            % r.t. ground truth
moving   =  start;
subpixel =  0.1;
[start_reg, output] =   efficient_imregister(fixed, moving, subpixel);
start4convmat       =   betterEdgeTaper(start_reg,option);                      % edge taper initial guess of g.t.
pncg_dI4convmat     =   start4convmat;
cg_dI4convmat       =   start4convmat;
%%%%%%%%%%%%%
% -------- ground truth comparison figure - start --------
drawComparisonFig(natureI,start, 0 ,'start','Nature',figPath);

% nature estimating
iterN =  1;     % one step for estimating nature
tolN  =  -inf;
HN    =  hessianMatrix(eye(imagesize));
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
        display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
end
        
for j = 1 : loopCompPerFrame
for i = 1 : methodsAmount
    method = option.method{i};
switch method
    case 'pncg'
    %%
%         startK = zeros(fsize);                   % initial guess of kernel, flat image with all NULL
        f_pncg = figure(114); clf(f_pncg); set(f_pncg,'visible','on')
        for jj = 1 : loopCompAllFrame            
            switch option.mode
                case 'compPerFrame'
                    k = j;
                case 'compAllFrame'
                    k = jj;
                otherwise
                    display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
            end
            % ##### special setup #####
            frame       =   clip(multiFrame{k},inf,0);
            natureK     =   multiKernel{k};                 % for error calculation
            
            % ##### estimate kernel #####
            %{%
            % comment out for non-blind deconvolution of ground truth, known PSFs, and learning theirs inverse
            % !!!!!!!! edge taper !!!!!!!!
            frameEdgeTaperred  =   betterEdgeTaper(frame,  option);
            switch option.blind
                case 'nb' % non blind
                    clear X
                    natureI4convmat   =   betterEdgeTaper(natureI,option);
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
                    
                case 'b'  % blind
                    clear X
                    switch option.img
                        case 'normal'
                            % !!!!!!!! normal image !!!!!!!!
                            X                   =   conv2MatOp(im2double(pncg_dI4convmat),fsize,shape);              % initial guess of convMtx X
                            startK              =   X'* frame ;                                                    % initial guess of kernel, b
                            frame4estiKernel    =   frameEdgeTaperred;
                        case 'gradient'
                            % !!!!!!!! gradient image !!!!!!!!
                            X                           =   conv2gradMat(im2double(pncg_dI4convmat),fsize,shape);    % initial guess of gradConvMtx X
                            frameGrad                   =   cell(1,2);
                            [frameGrad{1},frameGrad{2}] =   gradient(frameEdgeTaperred);
                            startK                      =   X'* frameGrad;                                         % initial guess of kernel, b
                            frame4estiKernel            =   frameGrad;
                    end
            end
            %
            option.plotFlag     =   1;
%             startK              =   startK./sum(vec(startK));
            startK              =   zeros(size(startK));
            [pncg_kernel, HK, errs_pncgK, ~, rerrs_pncgK] = deconv_pncg(X, frame4estiKernel, natureK, HK, iterK, startK, tolK, eta, option); % pncg
            pncg_kernel         =   preserveNorm(pncg_kernel);            % preserve energy norm of PSF
            % ----------- figure all S's -----------
            imgSize = size(HK.H);
            imgCells = cellImg(HK.s,imgSize);
            imgCelly = cellImg(HK.y,imgSize);
            tightSubplot(imgCells, [0,0], 'S', figPath, k)
            tightSubplot(imgCelly, [0,0], 'Y', figPath, k)
            % ----------- END all S's -----------
            % !!!!!!!! MEMLIM !!!!!!!!
            %{%
%             keyboard
            MEMLIM           =  option.MEMLIM;% size(HK.s,2);
            lambda           =  option.memoryStrength;
            [S,Y,Delta,GInv] =  purify(HK.s,HK.y,HK.delta,HK.Ginv0,MEMLIM,lambda);
            clear HK
            HK               =  hessianMatrix(eye(fsize)*scaler, S, Y, Delta, GInv, size(S,2)+1);
%             max(max(S'*Y - diag(1./diag(HK.Ginv0))))
            %}
            % ------- END MEMLIM -------
            % ----------- figure all Stilde's Ytilde's-----------
            imgSize = size(HK.H);
            imgCells = cellImg(HK.s,imgSize);
            imgCelly = cellImg(HK.y,imgSize);
            tightSubplot(imgCells, [0,0], 'Stilde', figPath, k)
            tightSubplot(imgCelly, [0,0], 'Ytilde', figPath, k)
            % ----------- END all S's -----------
            % !!!!!!!! DISCARD old Observations !!!!!!!!
            %{
            [S,Y,Delta]      =  discardObs(HK.s,HK.y,HK.delta,cutLine);
            clear HK
            HK               =  hessianMatrix(eye(fsize)*scaler, S, Y, Delta, nan, size(S,2)+1);
            cutLine          =  size(S,2)+1;
            %}
            % ------- END DISCARD -------
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,pncg_kernel,k,'pncg','Kernel',figPath);
            %}
            % ##### estimate nature #####
            if isfield(option,'blind') && strcmp(option.blind ,'b')
            % comment out for non-blind deconvolution of kernel, known ground truth, and learning its inverse
            option.plotFlag = 0;
            
            % !!!!!!!! blind !!!!!!!!
            clear Kpncg    
            Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
            clear HN
            HN = hessianMatrix(eye(imagesize));                          % normal CG step for g.t. estimation
            % --------- pncg step ---------
            if k == 1
                [pncg_dI, HN, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, start, tolN, eta, option); % pncg
            else
                [pncg_dI, HN, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, pncg_dI, tolN, eta, option); % pncg
            end
            % --------- gaussian step ---------
%             if k == 1
%                 [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,start,eta,option); % gaussian
%             else
%                 [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,pncg_dI,eta,option); % gaussian
%             end
            %
            % !!!!!!!! non blind with MEMLIM!!!!!!!!
% % %             option.plotFlag = 1;
% % %             clear Kpncg
% % %             Kpncg = conv2MatOp(im2double(natureK),imagesize,shape);  % convMtx of kernel, pncg
% % %             [pncg_dI, HN, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, 10, start, tolN, eta, option); % pncg
% % %             % ##### MEMLIM #####
% % %             MEMLIM = option.MEMLIM;% size(HK.s,2);
% % %             lambda = option.memoryStrength;
% % %             [S,Y,Delta,GInv] = purify(HN.s,HN.y,HN.delta,HN.Ginv0,MEMLIM,lambda);
% % %             clear HN
% % %             HN = hessianMatrix(eye(fsize)*scaler, S, Y, Delta, GInv, size(S,2)+1);
%             max(max(S'*Y - diag(1./diag(HK.Ginv0))))
            % !!!!!!!! END non blind with MEMLIM!!!!!!!!
            pncg_dI = clip(pncg_dI,inf,0);
            %%%%%%%%%%%%%
            fixed = natureI;                            % r.t. ground truth
            moving = pncg_dI;
            subpixel = 0.1;
            [pncg_dI, output] =   efficient_imregister(fixed, moving, subpixel);
            %%%%%%%%%%%%%
            pncg_dI4convmat =   betterEdgeTaper(pncg_dI,option);                      % edge taper every guess of g.t.
            % -------- ground truth comparison figure --------
            drawComparisonFig(natureI,pncg_dI,k,'pncg','Nature',figPath);
            
            % statitics 
            % all frame errors 
            errs_allframes_pncgK     = [errs_allframes_pncgK,  errs_pncgK];
            rerrs_allframes_pncgK    = [rerrs_allframes_pncgK, rerrs_pncgK];
            if k == 1
                errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(1), errs_pncgN(end)]; % residual, cg 
                rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(1),rerrs_pncgN(end)]; % relative error, cg
            else                
                errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(end)]; % residual, cg 
                rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(end)]; % relative error, cg
            end

            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, k, 'pncg', dre, figPath, 'debug', f_pncg)
            end
        end
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, numFrame, 'pncg', dre, figPath, 'latex');
            % -------- kernel frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_pncgK, rerrs_allframes_pncgK, numFrame, numFrame, 'pncg', dre, figPath, 'latex', f_pncg, 1);
            
            if exist('pncg_dI','var')
                I = pncg_dI;
            else
                I = 'non-blind';
            end
    case 'cg'
    %%
%         startK = zeros(fsize);                   % initial guess of kernel, flat image with all NULL
        f_cg = figure(115); clf(f_cg); set(f_cg,'visible','on')
        for jj = 1 : loopCompAllFrame            
            switch option.mode
                case 'compPerFrame'
                    k = j;
                case 'compAllFrame'
                    k= jj;
                otherwise
                    display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
            end
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            %{%
            % comment out for non-blind deconvolution of ground truth, known PSFs, and learning their inverse
            % !!!!!!!! edge taper !!!!!!!!
            frameEdgeTaperred  =   betterEdgeTaper(frame,  option);
            switch option.blind
                case 'nb' % non blind
                    clear X
                    natureI4convmat   =   betterEdgeTaper(natureI,option);
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
                    
                case 'b'  % blind
                    clear X
                    switch option.img
                        case 'normal'
                            % !!!!!!!! normal image !!!!!!!!
                            X                   =   conv2MatOp(im2double(cg_dI4convmat),fsize,shape);              % initial guess of convMtx X
                            startK              =   X'* frame ;                                                    % initial guess of kernel, b
                            frame4estiKernel    =   frameEdgeTaperred;
                        case 'gradient'
                            % !!!!!!!! gradient image !!!!!!!!
                            X                           =   conv2gradMat(im2double(cg_dI4convmat),fsize,shape);    % initial guess of gradConvMtx X
                            frameGrad                   =   cell(1,2);
                            [frameGrad{1},frameGrad{2}] =   gradient(frameEdgeTaperred);
                            startK                      =   X'* frameGrad;                                         % initial guess of kernel, b
                            frame4estiKernel            =   frameGrad;
                    end
            end
                    
            option.plotFlag             =   1;
            startK                      =   zeros(size(startK));
            [cg_kernel,errs_cg_kernel]  =   deconv_cg(X, frame4estiKernel, natureK, iterK, startK, tolK, eta, option); % cg
            cg_kernel                   =   preserveNorm(cg_kernel);                                                   % preserve energy norm of PSF
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,cg_kernel,k,'cg','Kernel',figPath);
            %}
%             keyboard
            % ##### estimate nature #####
            if isfield(option,'blind') && strcmp(option.blind ,'b')
            % comment out for non-blind deconvolution of kernel, known ground truth, and learning its inverse
            clear Kcg
            Kcg = conv2MatOp(im2double(cg_kernel),imagesize,shape);  % convMtx of kernel, cg
            option.plotFlag = 0;
            if k == 1
                [cg_dI, errs_cgN, ~, rerrs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, start, tolN, eta, option); % cg
            else
                [cg_dI, errs_cgN, ~, rerrs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, cg_dI, tolN, eta, option); % cg
            end
            % !!!!!!!! non blind !!!!!!!!
% % %             option.plotFlag = 1;
% % %             clear Kcg
% % %             Kcg = conv2MatOp(im2double(natureK),imagesize,shape);  % convMtx of kernel, cg
% % %             [cg_dI, errs_cgN, ~, rerrs_cgN] = deconv_cg(Kcg, frame, natureI, 10, start, tolN, eta, option); % cg
            % !!!!!!!! END non blind!!!!!!!!
            cg_dI = clip(cg_dI,1,0);
            %%%%%%%%%%%%%
            fixed = natureI;                            % r.t. ground truth
            moving = cg_dI;
            subpixel = 0.1;
            [cg_dI, output] =   efficient_imregister(fixed, moving, subpixel);
            %%%%%%%%%%%%%
            cg_dI4convmat   =   betterEdgeTaper(cg_dI,option);                        % edge taper every guess of g.t.
            %--------------- TOOL ---------------%
            % sum of cg_dI(it is more an image with sharp edges) and frame
            % only for illustration
            %
%             cg_dI_show = bsxfun(@times,cg_dI,(cg_dI>1e-5));
%             cg_dI_show = bsxfun(@times, cg_dI_show, sum(vec(frame))/sum(vec(cg_dI_show)));
%             cg_dI_show = cg_dI_show + frame;
%             cg_dI_show = bsxfun(@times, cg_dI_show, sum(vec(frame))/sum(vec(cg_dI_show)));
            %--------------- END  ---------------%
            % -------- ground truth comparison figure --------
            drawComparisonFig(natureI,cg_dI,k,'cg','Nature',figPath);

            % statitics 
            % all frame errors 
            if k == 1
                errs_allframes_cg  = [errs_allframes_cg, errs_cgN(1), errs_cgN(end)]; % residual, cg 
                rerrs_allframes_cg = [rerrs_allframes_cg,rerrs_cgN(1),rerrs_cgN(end)]; % relative error, cg
            else                
                errs_allframes_cg  = [errs_allframes_cg, errs_cgN(end)]; % residual, cg 
                rerrs_allframes_cg = [rerrs_allframes_cg,rerrs_cgN(end)]; % relative error, cg
            end
            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            drawAllFrameErrorFig(errs_allframes_cg, rerrs_allframes_cg, numFrame, k, 'cg', mpg, figPath, 'debug', f_cg);
            end
        end 
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_cg, rerrs_allframes_cg, numFrame, numFrame, 'cg', mpg, figPath, 'latex');
            %
            if exist('cg_dI','var')
                I = cg_dI;
            else
                I = 'non-blind';
            end
    case 'rl'
    %%
%         startK = ones(fsize)./prod(fsize);              % initial guess of kernel, flat image with uniform entries
        f_rl = figure(112); clf(f_rl); set(f_rl,'visible','on')
        for jj = 1 : loopCompAllFrame            
            switch option.mode
                case 'compPerFrame'
                    k = j;
                case 'compAllFrame'
                    k= jj;
                otherwise
                    display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
            end
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            % !!!!!!!! edge taper !!!!!!!!
            frame4estiKernel  =   betterEdgeTaper(frame,  option);
            % !!!!!!!! non blind !!!!!!!!
%             clear X
%             natureI4convmat   =   betterEdgeTaper(natureI,option);
%             X = conv2MatOp(im2double(natureI4convmat),fsize,shape);   % convMtx X of nature --- mfd
            %
            option.plotFlag = 1;
            [rl_kernel,errs_rl_kernel] = deconv_rl(X,frame4estiKernel,iterK,natureK,startK,eta,option); % richardson-lucy
            rl_kernel = preserveNorm(rl_kernel);            % preserve energy norm of PSF
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,rl_kernel,k,'rl','Kernel',figPath);
            
            % ##### estimate nature #####
            %{
            clear Krl
            Krl = conv2MatOp(im2double(rl_kernel),imagesize,shape);  % convMtx of kernel, gaussian
            option.plotFlag = 0;
            if k == 1
                [rl_dI,errs_rlN,rerrs_rlN] = deconv_rl(Krl,frame,iterN,natureI,start,eta,option); % gaussian
            else
                [rl_dI,errs_rlN,rerrs_rlN] = deconv_rl(Krl,frame,iterN,natureI,rl_dI,eta,option); % gaussian
            end
            clear X
            X = conv2MatOp(im2double(rl_dI),fsize,shape);   % new guess of convMtx X
            % -------- ground truth comparison figure --------
            drawComparisonFig(natureI,rl_dI,k,'rl','Nature',figPath);
            
            % statitics 
            % all frame errors
            errs_allframes_rl = [errs_allframes_rl,errs_rlN(end)]; % gaussian 
            rerrs_allframes_rl = [rerrs_allframes_rl,rerrs_rlN(end)]; % relative error, gaussian

            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            drawAllFrameErrorFig(errs_allframes_rl, rerrs_allframes_rl, numFrame, k, 'rl', ora, figPath, 'debug', f_rl);
            %}
        end            
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_rl, rerrs_allframes_rl, numFrame, numFrame, 'rl', ora, figPath, 'latex');
            %
            if exist('rl_dI','var')
                I = rl_dI;
            else
                I = 'non-blind';
            end
    case 'gaussian'
    %%
%         startK = ones(fsize)./prod(fsize);              % initial guess of kernel, flat image with uniform entries
        f_gaussian = figure(112); clf(f_gaussian); set(f_gaussian,'visible','on')
        for jj = 1 : loopCompAllFrame            
            switch option.mode
                case 'compPerFrame'
                    k = j;
                case 'compAllFrame'
                    k= jj;
                otherwise
                    display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
            end
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            % !!!!!!!! edge taper !!!!!!!!
            frame4estiKernel  =   betterEdgeTaper(frame,  option);
            % !!!!!!!! non blind !!!!!!!!
%             clear X
%             natureI4convmat   =   betterEdgeTaper(natureI,option);
%             X = conv2MatOp(im2double(natureI4convmat),fsize,shape);   % convMtx X of nature --- mfd
            %
            option.plotFlag = 1;
            [gaussian_kernel,errs_gaussian_kernel] = deconv_gaussian(X,frame4estiKernel,iterK,natureK,startK,eta,option); % gaussian
            gaussian_kernel = preserveNorm(gaussian_kernel);            % preserve energy norm of PSF
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,gaussian_kernel,k,'gaussian','Kernel',figPath);
            % ##### estimate nature #####
            %{
            clear Kgaussian
            Kgaussian = conv2MatOp(im2double(gaussian_kernel),imagesize,shape);  % convMtx of kernel, gaussian
            option.plotFlag = 0;
            if k == 1
                [gaussian_dI,errs_gaussianN,rerrs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,start,eta,option); % gaussian
            else
                [gaussian_dI,errs_gaussianN,rerrs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,gaussian_dI,eta,option); % gaussian
            end
            clear X
            X = conv2MatOp(im2double(gaussian_dI),fsize,shape);   % new guess of convMtx X
            % -------- ground truth comparison figure --------
            drawComparisonFig(natureI,gaussian_dI,k,'gaussian','Nature',figPath);
            gaussianNatureImg = figure;  set(gaussianNatureImg,'visible','off'),
            
            % statitics 
            % all frame errors
            errs_allframes_gaussian = [errs_allframes_gaussian,errs_gaussianN(end)]; % gaussian 
            rerrs_allframes_gaussian = [rerrs_allframes_gaussian,rerrs_gaussianN(end)]; % relative error, gaussian

            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            drawAllFrameErrorFig(errs_allframes_gaussian, rerrs_allframes_gaussian, numFrame, k,'gaussian', blu, figPath, 'debug', f_gaussian);
            %}
        end            
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_gaussian, rerrs_allframes_gaussian, numFrame, numFrame,'gaussian', blu, figPath, 'latex');
            %
            if exist('gaussian_dI','var')
                I = gaussian_dI;
            else
                I = 'non-blind';
            end
    otherwise
        display '[mbd.m]: option.method can be one of those: 'pncg', 'cg', 'rl', 'gaussian''
end
end
    

%##### set, label and save figure #####
saveResultFigure;

end
end
