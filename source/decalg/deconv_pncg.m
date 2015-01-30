function [pncg_dI,H,errs,tDeconv] = deconv_pncg(F,im,nature,H,iter,start,tol,eta,option)
% probabilistic solver, conjugate gradient

startup;
a = 1e-30; 

if nargin < 4
    H = hessianMatrix(eye(size(F'*im)));
end
if nargin < 5
    iter = 100;
end

if nargin < 6
    start = F'*im;
%     start = F'*(F*(F'*im));
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
x = start;
r0 = (F'*(F*x) +eta*((L*x)) + a*x) - b;
r = r0;
switch option.version
    case 'FH' % full H 
        if isempty(H.s) %H == eye(size(A)) %
            p = -r;
        else
            p = H*r;
        end
    case 'CG'
        p = -r;
    case 'SU' % speed up
        if isempty(H.s) %H == eye(size(A)) %
            p = -r;
        else
            p = H*r;
        end
    otherwise
        error('check option.version!')
end

epsl = 1e-30;
errs = nan(1,iter);
rerrs = nan(1,iter);
x = vec(x);
r = vec(r);
p = vec(p);
time = 1e-2;
residual = r;

%##### iteration #####
for k = 1:(iter + 1)  %numel(im)
    pncg_dI = reshape(x,imageSize);  
    im_residual = (F * pncg_dI - im); % 
    if norm(im_residual )==0
        errs(k) = 1e-20;
        rerrs(k) = 1e-20;
    else
        errs(k) = norm(im_residual)/ numel(im_residual); % average, absolute residual
        rerrs(k) = norm(pncg_dI - nature)/ norm(nature); % relative error ||x - hat(x)|| / ||x||
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
    
    if norm(im_residual) < numel(im_residual)*tol
        disp('==> solution found!')
        break
    else
        if k == (iter + 1)
            pncg_dI = reshape(x,imageSize);
            break
        end
        tStart = tic;
        % ###################
        q = vec((F'*(F*(reshape(p,imageSize)))  + eta*(L*reshape(p,imageSize)) + a*reshape(p,imageSize))); % A*p register
%         q = vec((F'*(F*(reshape(p,imageSize))) + a*reshape(p,imageSize))); % A*p register, without Tikhonov
        alpha = -(p'*q + epsl)\(p'*r);
        
        s = alpha*p;           % s_i <-- x_i+1 - x_i
        y = alpha*q;           % y_i <-- A*s_i
        
        switch option.version
            case 'FH'
%                 Hy = vec(H*reshape(y,imageSize));
%                 delta = s - Hy;        % delta_i <-- s_i - H_i*y_i
                delta = s - H.scale*y;          % ########29.01#########
            case 'CG'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'SU'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            otherwise
                error('check option.version!')
        end
             
        x = x + s;              % x_i+1 <-- x_i - alpha*p_i
        r = r + y;              % r_i+1 <-- r_i - A*alpa*p_i
        
%         den = s'*y;
%         H = H + (den + epsl)\(delta*s') + (den + epsl)\(s*delta') - ((den)^2 + epsl)\(s*delta'*y*s');% H_i+1 <-- H_i + (update)
        H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
        
        p_1 = p;
        switch option.version
            case 'FH'
                if k == 1
                    p = vec(H*(reshape(r,imageSize)));  % p <-- H*(A*x-b) = H*r;
                else
                    p = vec(H*(reshape(r,imageSize)));  % p <-- H*(A*x-b) = H*r;
%                     p = vec(Hy + p + H.*r);
                end
            case 'CG'
                p = r + H.*r;
            case 'SU'
                p = r + H.*r;
            otherwise
                error('check option.version!')
        end
        
        % ###################
        tElapsed = toc(tStart);
        % ########29.01###########
%         keyboard
%         buildH;
% orthogonal
        residual = [residual,r];
        display(sprintf('orth residual: %d', residual(:,end-1)'*residual(:,end)));
% conjugate
        conj = p_1'*vec(vec((F'*(F*reshape(p,imageSize))))  + vec(eta*(L*reshape(p,imageSize))) + a*p)
        time = [time;time(end)+tElapsed];
    end
    
end
tDeconv = time(end);
        % ########29.01###########
        buildH;
        [U,D] = eig(H_mtx);
        figure(6), imagesc(log10(abs(real(U'*U)))), colormap gray, axis image
        figure(7), imagesc(log10(diag(sort(diag(abs(real(D))),'ascend')))),    colormap gray, axis image
        print -depsc2 eigenspectrum_img
        keyboard
pncg_dI = clip(pncg_dI,1,0);


%##### figure #####
%----- main curves -----
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
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
hData = plot(errs, 'Color',dre); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
f13=figure(13); set(f13,'visible','off');
hData = plot(rerrs,'Color',dre); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
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