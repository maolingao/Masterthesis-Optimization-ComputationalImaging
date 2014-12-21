function [pncg_dI,H,errs,tDeconv] = deconv_pncg(F,im,nature,H,iter,start,tol,eta,option)
% probabilistic solver, conjugate gradient

startup;
a = 1e-12; 

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
x=vec(x);
r = vec(r);
p = vec(p);
time = 1e-2;


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
    loglog(errs,'Color',dre)
    ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex')
    xlabel('$\#steps$','Interpreter','Latex')
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
        alpha = (p'*q + epsl)\(p'*r);
        
        s = -alpha*p;           % s_i <-- x_i+1 - x_i
        y = -alpha*q;           % y_i <-- A*s_i
        
        switch option.version
            case 'FH'
                Hy = vec(H*reshape(y,imageSize));
                delta = s - Hy;        % delta_i <-- s_i - H_i*y_i
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
                    p = vec(Hy + p + H.*r);
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
% orthogonal
        orth = p_1'*vec(vec((F'*(F*reshape(p,imageSize))))  + vec(eta*(L*reshape(p,imageSize))) + a*p)
        time = [time;time(end)+tElapsed];
    end
    
end
tDeconv = time(end);
pncg_dI = clip(pncg_dI,1,0);


%##### figure #####
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
%----- main curves -----
% for debug
fclk = figure(14);  set(fclk,'visible','on'),subplot(121),loglog(time,errs,'Color',dre),hold on, 
subplot(122), set(fclk,'visible','on'),loglog(time,rerrs,'Color',dre),hold on
fstp = figure(15); set(fstp,'visible','on'),subplot(121),loglog(errs,'Color',dre),hold on, 
subplot(122), set(fstp,'visible','on'),loglog(rerrs,'Color',dre),hold on
% for latex
f10=figure(10); set(f10,'visible','off'),loglog(time,errs,'Color',dre),hold on, 
f12=figure(12); set(f12,'visible','off'),loglog(time,rerrs,'Color',dre),hold on
f11=figure(11); set(f11,'visible','off'),loglog(errs,'Color',dre),hold on, 
f13=figure(13); set(f13,'visible','off'),loglog(rerrs,'Color',dre),hold on
%
%----- image evolution and residual curve -----
figPath = '/is/ei/mgao/Documents/thesis/notes/meetingReport/fig';
%
f4 = figure(4); set(f4,'visible','on')
filename = ['deconv_pncg_with_curve','.tikz'];
filename = fullfile(figPath,filename);
  matlab2tikz(filename,'standalone',true,...
      'width','\fwidth','height','\fheight',...
      'parseStrings',false,...
      'extraAxisOptions',...
    {'xlabel near ticks','ylabel near ticks','scale only axis',...
    'label style={font=\footnotesize}', ...
    'legend style={font=\tiny}'...
    'title style={font=\footnotesize}'...
    'xticklabel style={font=\footnotesize}','yticklabel style={font=\footnotesize}'},...
    'showInfo', false);
%----- pncg deconved image -----
f_pncg = figure; set(f_pncg,'visible','off');
imagesc(clip(pncg_dI,1,0)); axis equal, colormap(gray)
title('my pncg')
filename = ['deconv_pncg','.tikz'];
filename = fullfile(figPath,filename);
  matlab2tikz(filename,'standalone',true,...
      'width','\fwidth','height','\fheight',...
      'parseStrings',false,...
      'extraAxisOptions',...
      {'title style={font=\small}','hide x axis', 'hide y axis'},...
      'showInfo', false);
end