function I = mfd_cu_gt(multiFrame, F, start, iter, natureI, multiKernel, option)
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
errs_allframes_pncgK = []; rerrs_allframes_pncgK = [];
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
scaler =  1e0;
H      =  hessianMatrix(eye(fsize)*scaler);
time  = 0;
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
            
            % !!!!!!!! normal image !!!!!!!!
            X                 =   conv2MatOp(im2double(natureK),imagesize,shape);   
            startK            =   X'* frame ;                                          
                    
            option.plotFlag     =   1;
            
            [pncg_gt, H, data_pncgK] = deconv_pncg(X, frame, natureI, H, iter, startK, tol, eta, option); % pncg
            pncg_gt         =   preserveNorm(pncg_gt);            % preserve energy norm of PSF
            % ###### psf residual curve ######
            data.errs(:,k) = data_pncgK.errs;
            data.rerrs(:,k) = data_pncgK.rerrs;
            data.time(:,k) = data_pncgK.time;
            figure(200), set(gcf,'visible','on')
            subplot(121), loglog(1:length(data_pncgK.errs),data_pncgK.errs); hold on
            subplot(122), loglog(data_pncgK.time,data_pncgK.errs); hold on
            
            % !!!!!!!! MEMLIM !!!!!!!!
            MEMLIM           =  option.MEMLIM;
            lambda           =  option.MEMSTR;
            % ################################# %
            [S,Y,Delta,GInv] =  purify(H.s,H.y,H.delta,MEMLIM,lambda);
            clear H
            H                =  hessianMatrix(eye(fsize)*scaler, S, Y, Delta, [] , []);
            % ################################# %
            % -------- kernel comparison figure --------
            drawComparisonFig(natureK,pncg_gt,k,'pncg','cmp',figPath);
       
            % statitics 
            % all frame errors - kernel
            time                     = timeLabel(time, data_pncgK.time);
            errs_allframes_pncgK     = [errs_allframes_pncgK,  data_pncgK.errs];
            rerrs_allframes_pncgK    = [rerrs_allframes_pncgK, data_pncgK.rerrs];
        end
            % -------- ground truth frame error figure --------
            % for latex
            % residual error & relative error
%             drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, numFrame, 'pncg', dre, figPath, 'latex');
            % -------- kernel frame error figure --------
            % for latex
            % residual error & relative error
            drawAllFrameErrorFigKernel(errs_allframes_pncgK, rerrs_allframes_pncgK, time, numFrame, numFrame, 'pncg', dre, figPath,  'latex')
            if exist('pncg_dI','var')
                I = pncg_dI;
            else
                I = pncg_gt;
            end
            
end 
switch option.version
    case 'CG'
        save('data_nb_cg.mat','data');
    case 'FH'
        save('data_nb_pcg.mat','data');
end
%##### set, label and save figure #####
saveResultFigure;

end
end
