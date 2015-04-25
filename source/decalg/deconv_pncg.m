function [pncg_dI, H, data] = deconv_pncg(F, im, nature, H, iter, start, tol, eta, option)
% probabilistic solver, conjugate gradient
% written for mfd_test script, to verify idea of multi-frame deconv(update of H with class hessianMatrix)

% example - input start = x, output must stay at start
startup;
a = 1e-6; 

if nargin < 4
    H = hessianMatrix(eye(size(F'*im)));
end
if nargin < 5
    iter = 100;
end

if nargin < 6
    start = F'*im;
end
if nargin < 7
    tol = 10^-6;
end
if nargin < 8
    eta = 1e-2;
end
if nargin < 9
    option.version = 'FH';
    option.figPath = '/is/ei/mgao/figure2drag';
    option.color = 'dre';
    option.LineStyle = '-';
end
if ~isfield(option,'color')
    color = dre;
else
    color = eval(option.color);
end
if ~isfield(option,'LineStyle')
    linestyle = '-';
else
    linestyle = option.LineStyle;
end
if ~isfield(option,'LineWidth')
    linewidth = 2;
else
    linewidth = option.LineWidth;
end

b           =   F'*im;
imageSize   =   size(b);
start       =   clip(start,inf,0);
%##### Tikhonov #####
l = [0 -1 0
     -1 4 -1
     0 -1 0]; % laplacian matrix WORK! 
L = conv2MatOp(l,imageSize,'same');
% eta = 0.01;

% ##### setup #####
%{%
switch option.version
    case 'FH'
%         M       =   hessianMatrix(H.H,H.s,H.y,H.delta); % preconditioner
        M       =   hessianMatrix(H.H,H.s,H.y,H.delta,H.R,H.D); % preconditioner
        x       =   start;
        r0      =   (F'*(F*x) + eta*((L*x)) + a*x) - b;
        r       =   r0;
%         p       =   M*r0;
        assert(isprop(H,'s') && isprop(H,'R'));
        if isempty(H.s) && isempty(H.R)
            p = -vec(r);
        else
            p = M*r;
        end
        if ~isfield(option,'color')
            switch option.img
                case 'normal'
                    color   =   blu;
                case 'gradient'
                    color   =   dre;
            end
        end
    case 'CG'
        M   = hessianMatrix(H.H); % preconditioner
        x       =   start;
        r0      =   (F'*(F*x) +eta*((L*x)) + a*x) - b;
        r       =   r0;
        p       =   M*r0;
        switch option.img
            case 'normal'
                color   =   blu;
            case 'gradient'
                color   =   mpg;
        end
    otherwise
        error('check option.version in deconv_pncg.m')
end
%}        
epsl         =  1e-30;
errs         =  nan(1,iter);
rerrs        =  nan(1,iter);
x            =  vec(x);
r            =  vec(r);
time         =  1e-2;
residual     =  r;
errRelChange =  nan;

%##### iteration #####
for k = 1 : (iter + 1)  %numel(im)
    pncg_dI         =   reshape(x,imageSize);  
    % -----------------------------------------------
    % if solving f, regularize kernel f
%     kernelSize = min(F.xsize, F.fsize);
%     if unique(abs(kernelSize - size(pncg_dI)) > abs(max(F.xsize, F.fsize) - size(pncg_dI)))
%         NOP;                                    % current solving x, BOP    
%     else
%         pncg_dI       =   lowerBound(pncg_dI);          % current solving f, low bound f
% %         pncg_dI       =   preserveNorm(pncg_dI);        % preserve energy norm of f
%     end
    % -----------------------------------------------
    % residual error
    im_residual     =   betterMinus(F * pncg_dI, im); % 
    % -----------------------------------------------
    % crop away edges
    kernelSize      =   min(F.xsize, F.fsize);
    corpMarginSize  =   kernelSize;
    Pim             =   patimat('same',size(im_residual),corpMarginSize,0);
    im_residual     =   Pim'*im_residual;
    % -----------------------------------------------
    % register images r.t. ground truth 
    fixed           =   nature;                            % r.t. ground truth
    moving          =   pncg_dI;
    subpixel        =   1;
    [pncg_dI_reg, output] = efficient_imregister(fixed, moving, subpixel);
    % -----------------------------------------------
    % absolute error
    errorabso       =   pncg_dI_reg - nature;
    if unique(abs(kernelSize - size(pncg_dI)) > abs(max(F.xsize, F.fsize) - size(pncg_dI)))
        errorabso   =   Pim'*errorabso;        % current solving x, crop    
        natureCrop  =   Pim'*nature;
    else
        errorabso   =   errorabso;             % current solving f, NOT crop
        natureCrop  =   nature;
    end
    % -----------------------------------------------
    % average residual & relative error of ground truth
    if norm(im_residual,'fro' )==0
        errs(k)     =   1e-20;
        rerrs(k)    =   1e-20;
    else            
        errs(k)     =   (norm(im_residual,'fro') / numel(im_residual)); % average, absolute residual
        rerrs(k)    =   (norm(errorabso,'fro') / norm(natureCrop,'fro')); % relative error ||x - hat(x)|| / ||x||
    end
    % -----------------------------------------------
    % plot 
    f4 = figure(4); subplot(131)
    imagesc(clip(nature,1,0)); axis image,colormap(gray)
    title('ground truth')
    drawnow          
    subplot(132)
    imagesc(clip(pncg_dI,1,0)); axis image,colormap(gray)
    title(sprintf('my pncg - iteration %d/%d',k,iter + 1))
    drawnow          
    subplot(133)
    hData = loglog(errs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth);
    hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
    hXLabel = xlabel('$\#steps$', 'Interpreter','Latex');
    thisFigure;   
    drawnow 
    % -----------------------------------------------
    % stop creterien 1 : solution found
    if norm(im_residual,'fro') < numel(im_residual)*tol
        disp('==> solution found!')
        break
    else
    % -----------------------------------------------
    % stop creterien 2 : error decreasing tiny or even increasing
        if k > 3
            errRelChange    =   errs(2:k) - errs(1:k-1);
            errRelChange    =   sum(errRelChange(k-3:k-1)) / sum(errs(k-3:k)) * 4/3 ;
        end
    % -----------------------------------------------
    % stop creterien 3 : iteration number reached
        if k == (iter + 1) || errRelChange >   inf % -1e-3 %
            break
        end
    % -----------------------------------------------
    % main calculation
        p = bsxfun(@rdivide,p,sqrt(sum(p.^2)));     % normalize every direction, columnwise, this helps stabilize
        tStart = tic;
        % ###################
        q           =   nan(size(p));
        for i = 1 : size(p,2)
            q(:,i)  =   vec((F'*(F*(reshape(p(:,i),imageSize)))  + eta*(L*reshape(p(:,i),imageSize)) + a*reshape(p(:,i),imageSize))); % A*p register
        end
%         alpha       =   - (p'*q + epsl)\(p'*r);
        alpha       =   - pinv(p'*q) * (p'*r);
        
        s           =   p*alpha;                % s_i <-- x_i+1 - x_i
        y           =   q*alpha;                % y_i <-- A*s_i
        s_length    =   norm(s);
        s           =   s ./ s_length;          % for numerical stability, up scale s and y
        y           =   y ./ s_length;
        
        switch option.version
            case 'FH'
                delta   =   s - y;              % delta_i <-- s_i - H_i*y_i
            case 'CG'
                delta   =   s - y;              % delta_i <-- s_i - H_i*y_i
            otherwise
                error('check option.version in deconv_pncg.m')
        end
             
        x           =   x + p*alpha; % x_i+1 <-- x_i - alpha*p_i  =  x + s;
        r           =   r + q*alpha; % r_i+1 <-- r_i - A*alpa*p_i  =  r + y;
       
        H           =   plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
        
        p_1         =   p;
        switch option.version
            case 'FH'
                p   =   vec(H*(reshape(r,imageSize)));  % p <-- H*(A*x-b) = H_i+1 * r_i+1
%                 g   =   pinv(H.s'*H.y);
%                 p   =   r - H.s * g * (H.y' * r); % this ensure conjugacy
            case 'CG'
                p   =   r + H.*r;
            otherwise
                error('check option.version in deconv_pncg.m')
        end
        % ###################
        tElapsed = toc(tStart);
        time     = [time, time(end)+tElapsed];
    % -----------------------------------------------
    % display orthogonality and conjugacy 
    % orthogonal      
        residual = [residual,r]; 
        display(sprintf('orth residual: %d', residual(:,end-1)'*residual(:,end)))
    % conjugate
        conj = p_1'*vec(vec((F'*(F*reshape(p,imageSize))))  + vec(eta*(L*reshape(p,imageSize))) + a*p)
    end
    
end
tDeconv = time(end);            % total time for current deconvolution task
pncg_dI = clip(pncg_dI,1,0);    



%##### figure #####
%----- main curves -----
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
% rerrs = rerrs./max(rerrs); % normalize rel. error
if isfield('option', 'plotFlag') && option.plotFlag == 1
    
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121), hData = loglog(time ,errs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); thisFigure; hold on
subplot(122), hData = loglog(time,rerrs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); thisFigure; hold on
fstp = figure(15); set(fstp,'visible','on'),
subplot(121), hData = loglog( errs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); thisFigure; hold on
subplot(122), hData = loglog(rerrs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); thisFigure; hold on
% for latex

    etaVec = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0];
    SNRVec = 10:10:60;
    switch option.mode
        case 'SNR'
            trgVec = SNRVec;
            for k = 1 : length(trgVec)
                content{k} = strcat('SNR', sprintf(' %gdB',trgVec(k)));
            end
        case 'eta'
            trgVec = etaVec;
            for k = 1 : length(trgVec)
                content{k} = strcat('$\eta$', sprintf(' %g',trgVec(k)));
            end
        otherwise
            content = [];
    end
% content{1} = 'gradient image';
% content{2} = 'normal image';

f10=figure(10); set(f10,'visible','off');
hData = loglog(time, errs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); 
%     hLegend = legend(content); 
%     set(hLegend,'Interpreter','Latex','Location','northeast');
    hXLabel = xlabel('$time(sec)$');
    hYLabel = ylabel('$residual\ error$');
axis tight; thisFigure; hold on
f12=figure(12); set(f12,'visible','off');
hData = loglog(time,rerrs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); 
%     hLegend = legend(content); 
%     set(hLegend,'Interpreter','Latex','Location','northeast');
    hXLabel = xlabel('$time(sec)$');
    hYLabel = ylabel('$relative\ error$');
axis tight; thisFigure; hold on
f11=figure(11); set(f11,'visible','off');
hData = loglog(errs, 'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); 
%     hLegend = legend(content); 
%     set(hLegend,'Interpreter','Latex','Location','northeast');
    hXLabel = xlabel('$\#steps$');
    hYLabel = ylabel('$residual\ error$');
% set(gca,'Yscale','log'), 
axis tight; thisFigure; hold on 
f13=figure(13); set(f13,'visible','off');
hData = loglog(rerrs,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth); 
%     hLegend = legend(content); 
%     set(hLegend,'Interpreter','Latex','Location','northeast');
    hXLabel = xlabel('$\#steps$');
    hYLabel = ylabel('$relative\ error$');
% set(gca,'Yscale','log'), 
axis tight; thisFigure; hold on 
end
%
%----- image evolution and residual curve -----
figPath = option.figPath;

f4 = figure(4); set(f4,'visible','on')
filename = 'deconv_pncg_with_curve';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

%----- pncg deconved image -----
f_pncg = figure; set(f_pncg,'visible','off');
imagesc(clip(pncg_dI,1,0)); axis image off, colormap(gray)
title('pncg')
filename = 'deconv_pncg';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf
%----- relative error -----
% figure(34), set(gcf,'visible','on');
% plot(rerrs,'Color',color,'LineStyle',linestyle), hold on;
% h = legend('$\eta\ 0$','$\eta\ 0.001$','$\eta\ 0.01$','$\eta\ 0.1$','$\eta\ 1$','Location','NorthWest');
% set(h,'Interpreter','latex')
% % legend('SNR 10','SNR 20','SNR 30','SNR 40','SNR 50','Location','NorthWest')
% ylabel('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex') 
% % ylabel('$Relative error$','Interpreter','Latex')
% xlabel('$\#steps$','Interpreter','Latex')
% filename = 'deconv_pncg_relativeError';
% filename = fullfile(figPath,filename);
% print(gcf, '-depsc2', filename)
data.errs = errs;
data.time = time;
data.rerrs = rerrs;
end