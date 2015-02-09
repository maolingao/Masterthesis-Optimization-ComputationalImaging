function [pncg_dI, H, errs, tDeconv, rerrs] = deconv_pncg(F, im, nature, H, iter, start, tol, eta, option)
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

b = F'*im;
imageSize = size(b);
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
%         keyboard
        M = hessianMatrix(H.H,H.s,H.y,H.delta,H.Ginv0,H.i); % preconditioner
        x = start;
        r0 = (F'*(F*x) +eta*((L*x)) + a*x) - b;
        r = r0;
%         keyboard
        p = M*r0;
%         p = -r0;
        color = dre;
    case 'CG'
        M = hessianMatrix(H.H); % preconditioner
%         M = hessianMatrix(H.H,H.s,H.y,H.delta,H.Ginv0,H.i); % preconditioner
        x = start;
        r0 = (F'*(F*x) +eta*((L*x)) + a*x) - b;
        r = r0;
        p = M*r0;
        color = mpg;
    otherwise
        error('check option.version in pncg_Hmfd')
end
%}        
epsl = 1e-30;
errs = nan(1,iter);
rerrs = nan(1,iter);
x = vec(x);
r = vec(r);
p = vec(p);
time = 1e-2;
residual = r;
errRelChange = nan;

%##### iteration #####
% keyboard
for k = 1 : (iter + 1)  %numel(im)
    pncg_dI = reshape(x,imageSize);  
    im_residual = (F * pncg_dI - im); % 
    % crop away edge
%     keyboard
    kernelSize = min(F.xsize, F.fsize);
    corpMarginSize = kernelSize;
    Pim = patimat('same',size(im_residual),corpMarginSize,0);
    im_residual = Pim'*im_residual;
%           im_residual = (cg_dI - double(nature)) ;
    % absolute error
    fixed = nature;                            % r.t. ground truth
    moving = pncg_dI;
    subpixel = 1;
    [pncg_dI_reg, output] = efficient_imregister(fixed, moving, subpixel);
    errorabso = pncg_dI_reg - nature;
    if unique(abs(kernelSize - size(pncg_dI)) > abs(max(F.xsize, F.fsize) - size(pncg_dI)))
%         keyboard
%         figure(5), imagesc(im_residual), colormap gray, axis image
        errorabso = Pim'*errorabso;        % current solving x, crop    
        natureCrop = Pim'*nature;
    else
        errorabso = errorabso;             % current solving f, NOT crop
        natureCrop = nature;
    end
    % average residual & relative error of ground truth
    if norm(im_residual,'fro' )==0
        errs(k) = 1e-20;
        rerrs(k) = 1e-20;
    else            
        errs(k) = (norm(im_residual,'fro') / numel(im_residual)); % average, absolute residual
        rerrs(k) = (norm(errorabso,'fro') / norm(natureCrop,'fro')); % relative error ||x - hat(x)|| / ||x||
    end

    f4 = figure(4);
    subplot(121)
    imagesc(clip(pncg_dI,1,0)); axis image,colormap(gray)
    title(sprintf('my pncg - iteration %d/%d',k,iter + 1))
    drawnow          

    subplot(122)
    hData = loglog(errs,'Color',color,'LineStyle',linestyle);
    hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
    hXLabel = xlabel('$\#steps$', 'Interpreter','Latex');
    thisFigure;   
    drawnow 
    %

    if norm(im_residual,'fro') < numel(im_residual)*tol
        disp('==> solution found!')
        break
    else
        % stop creterien ########
        if k > 3
%             keyboard
            errRelChange = errs(2:k) - errs(1:k-1);
            errRelChange = sum(errRelChange(k-3:k-1)) / sum(errs(k-3:k)) * 4/3 ;
        end
        
        if k == (iter + 1) || errRelChange > -1e-3
%             keyboard
            pncg_dI = clip(reshape(x,imageSize),1,0);
            break
        end
        tStart = tic;
        % ###################
        q = vec((F'*(F*(reshape(p,imageSize)))  + eta*(L*reshape(p,imageSize)) + a*reshape(p,imageSize))); % A*p register
        alpha = - (p'*q + epsl)\(p'*r);
        
        s = alpha*p;           % s_i <-- x_i+1 - x_i
        y = alpha*q;           % y_i <-- A*s_i
        
        switch option.version
            case 'FH'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'CG'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            otherwise
                error('check option.version in pncg_Hmfd')
        end
             
        x = x + s;              % x_i+1 <-- x_i - alpha*p_i
        r = r + y;              % r_i+1 <-- r_i - A*alpa*p_i        
       
        H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
% % %         if s'*y > -1e-15
% % %             if option.frame > inf
% % %             % eliminate the component in new tuple [s,y,delta] which are
% % %             % already contained in previous spanned Krylov subspace.
% % %                 [S,Y,Delta] = purify(H,s,y,delta); 
% % % %                 keyboard
% % %                 clear H
% % %                 H = hessianMatrix(M.H,S,Y,Delta,(M.i+1));
% % %             else
% % %                 H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
% % %             end
% % %         else
% % %                     keyboard
% % %             display 'direction changes too small.'
% % %             break
% % %         end
        
        p_1 = p;
        switch option.version
            case 'FH'
                p = vec(H*(reshape(r,imageSize)));  % p <-- H*(A*x-b) = H_i+1 * r_i+1
            case 'CG'
                p = r + H.*r;
            otherwise
                error('check option.version in pncg_Hmfd')
        end
        % ###################
        tElapsed = toc(tStart);
%         keyboard
%         buildH;
% orthogonal      
        residual = [residual,r]; 
        display(sprintf('orth residual: %d', residual(:,end-1)'*residual(:,end)))
% conjugate
        conj = p_1'*vec(vec((F'*(F*reshape(p,imageSize))))  + vec(eta*(L*reshape(p,imageSize))) + a*p)
        time = [time;time(end)+tElapsed];
    end
    
end
tDeconv = time(end);
% keyboard
pncg_dI = clip(pncg_dI,1,0);

%
%
%##### figure #####
%----- main curves -----
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
if option.plotFlag == 1
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121), hData = loglog(time ,errs,'Color',dre); thisFigure; hold on
subplot(122), hData = loglog(time,rerrs,'Color',dre); thisFigure; hold on
fstp = figure(15); set(fstp,'visible','on'),
subplot(121), hData = loglog( errs,'Color',dre); thisFigure; hold on
subplot(122), hData = loglog(rerrs,'Color',dre); thisFigure; hold on
% for latex
f10=figure(10); set(f10,'visible','off');
hData = loglog(time, errs,'Color',dre); 
axis tight; thisFigure; hold on
f12=figure(12); set(f12,'visible','off');
hData = loglog(time,rerrs,'Color',dre); 
axis tight; thisFigure; hold on
f11=figure(11); set(f11,'visible','off');
hData = loglog(errs, 'Color',dre); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
f13=figure(13); set(f13,'visible','off');
hData = loglog(rerrs,'Color',dre); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
end
%
%----- image evolution and residual curve -----
figPath = option.figPath;

f4 = figure(4); set(f4,'visible','on')
filename = 'deconv_pncg_with_curve';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% keyboard
%----- pncg deconved image -----
f_pncg = figure; set(f_pncg,'visible','off');
imagesc(clip(pncg_dI,1,0)); axis image off, colormap(gray)
title('pncg')
filename = 'deconv_pncg';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
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
end