function [I, err]= mbd(multiFrame, F, start, iterK, natureI, multiKernel, option)
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
if isfield('option','aux'); aux = option.aux; else aux = 'noCmp'; end

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

% ##### general setup #####
numFrame = numel(multiFrame);
for i = 1 : numFrame                            % register observation images
    fixed           =   natureI;                            % r.t. ground truth
    moving          =   multiFrame{i};
    subpixel        =   1;
    [multiFrame{i}, output] = efficient_imregister(fixed, moving, subpixel);
    %
    multiKernel{i}  =  center(multiKernel{i});   % center all kernel, only for error analysis
end
etax = option.etax;                               % ### <--- regularization parameter
etaf = option.etaf;                               % ### <--- regularization parameter
% kernel estimation
tolK                =   option.tolK;
scaler              =   1e0;
HK                  =   hessianMatrix(eye(fsize)*scaler);
start               =   clip(start,inf,0);
start(start<1e-7)   =   0;
timeK               =   0;
% start  =  start./max(vec(start)) * 1e-3;
%%%%%%%%%%%%%
fixed               =   natureI;                            % r.t. ground truth
moving              =   start;
subpixel            =   1;
[start_reg, output] =   efficient_imregister(fixed, moving, subpixel);
tau = 0.3;   m      =   mask(start_reg, tau);
start4convmat       =   betterEdgeTaper(start_reg, option);                      % edge taper initial guess of g.t.
pncg_dI4convmat     =   start4convmat;
cg_dI4convmat       =   start4convmat;
%%%%%%%%%%%%%
% -------- ground truth comparison figure - start --------
drawComparisonFig(natureI,start, 0 ,'start','Nature',figPath,aux);

% nature estimation
iterN               =   1;     % one step for estimating nature
tolN                =   -inf;
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
            
            % ################################# estimate kernel #################################
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
                    % ^^^^^^^^^^^
                    % rmap mask
                    if ~isfield(option,'flagMask') || option.flagMask == 0                       
                        m = ones(size(m));
                    end
%                     figure, set(gcf,'visible','off')
%                     imagesc(m), title('m'); colormap gray, axis image off
%                     filename = sprintf('mask_%d',k);
%                     filename = fullfile(figPath,filename);
%                     print(gcf, '-depsc2', filename)
%                     close gcf;
                    % ^^^^^^^^^^^
                    switch option.img
                        case 'normal'
                            % !!!!!!!! normal image !!!!!!!!
                            X                   =   conv2MatOp(im2double(pncg_dI4convmat),fsize,shape);              % initial guess of convMtx X
                            startK              =   X'* frame ;                                                      % initial guess of kernel, b
                            frame4estiKernel    =   frameEdgeTaperred;
                        case 'gradient'
                            % !!!!!!!! gradient image !!!!!!!!
                            X                           =   conv2gradMat(im2double(pncg_dI4convmat),fsize,shape,m);  % initial guess of gradConvMtx X
                            frameGrad                   =   cell(1,2);
                            [frameGrad{1},frameGrad{2}] =   gradient(frameEdgeTaperred);
                            frameGrad                   =   cellfun(@(x)m.*x, frameGrad, 'UniformOutput', false);
                            %^^^^^^^
                            % gradImg for deconvolve psf
%                             [gx,gy]   = gradient(pncg_dI4convmat);                  % gradient of x
%                             gx        = m .* gx;
%                             gy        = m .* gy;                
%                             figure, set(gcf,'visible','off');
%                             subplot(2,2,1), imagesc(gx),colormap gray, axis image off, title('X_u')
%                             subplot(2,2,3), imagesc(gy),colormap gray, axis image off, title('X_v')
%                             subplot(2,2,2), imagesc(frameGrad{1}),colormap gray, axis image off, title('y_u')
%                             subplot(2,2,4), imagesc(frameGrad{2}),colormap gray, axis image off, title('y_v')
%                             filename = sprintf('gradImg4deconvPSF_%d',k);
%                             filename = fullfile(figPath,filename);
%                             print(gcf, '-depsc2', filename)
%                             close gcf;
                            %^^^^^^^
                            startK                      =   X'* frameGrad;                                           % initial guess of kernel, b
                            frame4estiKernel            =   frameGrad;
                    end
            end
            %
            option.plotFlag     =   1;
            [pncg_kernel, HK, errs_pncgK, clkK, rerrs_pncgK] = deconv_pncg(X, frame4estiKernel, natureK, HK, iterK, startK, tolK, etaf, option); % pncg
            pncg_kernel         =   preserveNorm(pncg_kernel);            % preserve energy norm of PSF
            
            % !!!!!!!! MEMLIM !!!!!!!!
            %{%
            MEMLIM           =  option.MEMLIM;
            lambda           =  option.MEMSTR;
            alpha            =  option.EXPOSTR;
            % ################################# %
            [S,Y,Delta,GInv] =  purify(HK.s,HK.y,HK.delta,MEMLIM,lambda);
            clear HK
            HK               =  hessianMatrix(eye(fsize)*scaler, S, Y, Delta, [] , []);
            % ------- END MEMLIM -------
            
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,pncg_kernel,k,'pncg','Kernel',figPath,aux);
            %}
            % ################################# estimate nature #################################
            if isfield(option,'blind') && strcmp(option.blind ,'b')
                % comment out for non-blind deconvolution of kernel, known ground truth, and learning its inverse
                option.plotFlag = 0;

                % !!!!!!!! blind !!!!!!!!
                clear Kpncg    
                Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
                clear HN
                HN = hessianMatrix(eye(imagesize));                          % normal CG step for g.t. estimation
                % --------- pncg step ---------
%                 if k == 1
%                     [pncg_dI, HN, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, start, tolN, eta, option); % pncg
%                 else
%                     [pncg_dI, HN, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, pncg_dI, tolN, eta, option); % pncg
%                 end
                % --------- gaussian step ---------
                if k == 1
                    [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,start,etax,option); % gaussian
                else
                    [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,pncg_dI,etax,option); % gaussian
                end
                
                pncg_dI     =   clip(pncg_dI,inf,0);
                %%%%%%%%%%%%%
                fixed       =   natureI;                            % r.t. ground truth
                moving      =   pncg_dI;
                subpixel    =   0.1;
                [pncg_dI, output]   =   efficient_imregister(fixed, moving, subpixel);
                %%%%%%%%%%%%%
                %^^^^^^^^^^^^^
                m = mask(pncg_dI,tau);
                %^^^^^^^^^^^^^
                pncg_dI4convmat     =   betterEdgeTaper(pncg_dI,option);                      % edge taper every guess of g.t.
                % -------- ground truth comparison figure --------
                drawComparisonFig(frame,pncg_dI,k,'pncg','Nature',figPath,aux);
                % statitics 
                % all frame errors - ground truth
                if k == 1
                    errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(1), errs_pncgN(end)]; % residual, pncg 
                    rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(1),rerrs_pncgN(end)]; % relative error, pncg
                else                
                    errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(end)]; % residual, pncg 
                    rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(end)]; % relative error, pncg
                end
                % -------- ground truth frame error figure --------
                % for debug
                drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, k, 'pncg', dre, figPath, 'debug', f_pncg)
            end            
            % statitics 
            % all frame errors - kernel
            timeK                    = timeLabel(timeK, clkK);
            errs_allframes_pncgK     = [errs_allframes_pncgK,  errs_pncgK];
            rerrs_allframes_pncgK    = [rerrs_allframes_pncgK, rerrs_pncgK];
        end
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, numFrame, 'pncg', dre, figPath, 'latex');
            % -------- kernel frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFigKernel(errs_allframes_pncgK, rerrs_allframes_pncgK, timeK, numFrame, numFrame, 'pncg', dre, figPath,  'latex')
            if exist('pncg_dI','var')
                I = pncg_dI;
            else
                I = 'non-blind';
            end
            err.res = errs_allframes_pncg;
            err.rel = rerrs_allframes_pncg;
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
            drawComparisonFig(natureK,cg_kernel,k,'cg','Kernel',figPath,aux);
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
            drawComparisonFig(frame,cg_dI,k,'cg','Nature',figPath,aux);

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
            err.res = errs_allframes_cg;
            err.rel = rerrs_allframes_cg;
    otherwise
        display '[mbd.m]: option.method can be one of those: 'pncg', 'cg''
end
end
    

%##### set, label and save figure #####
saveResultFigure;

end
end
