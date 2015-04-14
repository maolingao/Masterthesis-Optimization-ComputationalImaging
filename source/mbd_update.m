function [I, err] = mbd_update(multiFrame, F, start, iterK, natureI, multiKernel, option)
% mbd wrap up pncg

if nargin < 7
    option.figPath = '/is/ei/mgao/figure2drag';
    option.win = 'barthann';
    option.MEMLIM = 20;
end

startup;
figPath = option.figPath;
etax = option.etax;                               % ### <--- regularization parameter
etaf = option.etaf;                               % ### <--- regularization parameter

if isfield(option,'aux'); aux = option.aux; else aux = 'noCmp'; end

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

% ##### general setup #####
numFrame = numel(multiFrame);
for i = 1 : numFrame                            % register observation images
    fixed           =   natureI;                            % r.t. ground truth
    moving          =   multiFrame{i};
    subpixel        =   1;
    [multiFrame{i}, output] = efficient_imregister(fixed, moving, subpixel);
end

% kernel estimation
tolK                =   option.tolK;
scaler              =   1e0;
HK                  =   hessianMatrix(eye(fsize)*scaler);
start               =   clip(start,inf,0);
start(start<1e-7)   =   0;
timeK               =   0;

start4convmat       =   betterEdgeTaper(start, option);                      % edge taper initial guess of g.t.
pncg_dI4convmat     =   start4convmat;
% rmap mask
tau = 0.3;   m      =   mask(start, tau);
if ~isfield(option,'flagMask') || option.flagMask == 0                       
    m = ones(size(m));
end
% ground truth comparison figure - start
drawComparisonFig(natureI,start, 0 ,'start','Nature',figPath,aux);

% nature estimation
iterN               =   1;     % one step for estimating nature
tolN                =   -inf;

% iteration
%         startK = zeros(fsize);                   % initial guess of kernel, flat image with all NULL
f_pncg = figure(114); clf(f_pncg); set(f_pncg,'visible','on')

for k = 1 : numFrame
    % special setup
    frame       =   clip(multiFrame{k},inf,0);
    natureK     =   multiKernel{k};                 % for error calculation

    %% ^^^^^^^ PSF estimation ^^^^^^^
    % edge taper 
    frameEdgeTaperred  =   betterEdgeTaper(frame,  option);
    % gradient image
    X                           =   conv2gradMat(im2double(pncg_dI4convmat),fsize,shape,m);  % initial guess of gradConvMtx X
    frameGrad                   =   cell(1,2);
    [frameGrad{1},frameGrad{2}] =   gradient(frameEdgeTaperred);
    frameGrad                   =   cellfun(@(x)m.*x, frameGrad, 'UniformOutput', false);
    startK                      =   X'* frameGrad;                                           % initial guess of kernel, b
    frame4estiKernel            =   frameGrad;
    %
    option.plotFlag     =   1;
    [pncg_kernel, HK, errs_pncgK, clkK, rerrs_pncgK] = deconv_pncg(X, frame4estiKernel, natureK, HK, iterK, startK, tolK, etaf, option); % pncg
    pncg_kernel         =   preserveNorm(pncg_kernel);            % preserve energy norm of PSF
   
    % kernel comparison figure
    drawComparisonFig(natureK,pncg_kernel,k,'pncg','Kernel',figPath,aux);

    %% ^^^^^^^ MEMLIM ^^^^^^^
    MEMLIM           =  option.MEMLIM;
    lambda           =  option.MEMSTR;
    alpha            =  option.EXPOSTR;

    [S,Y,Delta,GInv] =  purify(HK.s,HK.y,HK.delta,MEMLIM,lambda);
    clear HK
    HK               =  hessianMatrix(eye(fsize)*scaler, S, Y, Delta, [] , []);

    %% ^^^^^^^ g.t. estimation ^^^^^^^
    option.plotFlag = 0;

    clear Kpncg    
    Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
    % --------- pncg step ---------
%     clear HN
%     HN = hessianMatrix(eye(imagesize));                          % normal CG step for g.t. estimation
%     if k == 1
%         [pncg_dI, ~, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, start, tolN, etax, option); % pncg
%     else
%         [pncg_dI, ~, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, pncg_dI, tolN, etax, option); % pncg
%     end
    % --------- gaussian step ---------
    if k == 1
        [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,start,etax,option); % gaussian
    else
        [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,pncg_dI,etax,option); % gaussian
    end

    pncg_dI             =   clip(pncg_dI,inf,0);
    pncg_dI4convmat     =   betterEdgeTaper(pncg_dI,option);                      % edge taper every guess of g.t.

    % ground truth comparison figure 
    drawComparisonFig(frame,pncg_dI,k,'pncg','Nature',figPath,aux);

    %% statitics 
    % all frame errors - ground truth
    if k == 1
        errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(1), errs_pncgN(end)]; % residual, pncg 
        rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(1),rerrs_pncgN(end)]; % relative error, pncg
    else                
        errs_allframes_pncg  = [errs_allframes_pncg, errs_pncgN(end)]; % residual, pncg 
        rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(end)]; % relative error, pncg
    end

    % ground truth frame error figure  
    % for debug
    drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, k, 'pncg', dre, figPath, 'debug', f_pncg)

    % statitics 
    % all frame errors - kernel
    timeK                    = timeLabel(timeK, clkK);
    errs_allframes_pncgK     = [errs_allframes_pncgK,  errs_pncgK];
    rerrs_allframes_pncgK    = [rerrs_allframes_pncgK, rerrs_pncgK];
end

% ground truth frame error figure 
% for latex
% residual error & relative error
drawAllFrameErrorFig(errs_allframes_pncg, rerrs_allframes_pncg, numFrame, numFrame, 'pncg', dre, figPath, 'latex');
% kernel frame error figure 
% for latex
% residual error & relative error
drawAllFrameErrorFigKernel(errs_allframes_pncgK, rerrs_allframes_pncgK, timeK, numFrame, numFrame, 'pncg', dre, figPath,  'latex')

% output
I = pncg_dI;

err.res = errs_allframes_pncg;
err.rel = rerrs_allframes_pncg;

%##### set, label and save figure #####
saveResultFigure;

