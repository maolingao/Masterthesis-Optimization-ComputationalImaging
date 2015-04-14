function [gaussian_dI,errs,rerrs] = deconv_gaussian(F,im,iter,nature,start,eta,option)
% Gaussian noise deconvolution
startup;

if nargin < 5
    start = F'*im;
end

if nargin < 6
    eta = 0;
end

if nargin < 7
    option.figPath = '/is/ei/mgao/figure2drag';
end

%##### setup #####
gaussian_dI = start; % custermized start guess; 
% gaussian_dI = nature; % custermized start guess; 

errs = nan(1,iter);
rerrs = nan(1,iter);
            
epsl = 1e-10;
time = 1e-2;
errRelChange = nan;
%##### Tikhonov #####
%{
l = [0 -1 0
     -1 4 -1
     0 -1 0]; % laplacian matrix WORK! 
lp = [0 0 0
     0 4 0
     0 0 0]; % positive laplacian matrix
ln = lp - l; % negative laplacian matrix
Lp = conv2MatOp(lp,imageSize,F.shape);
Ln = conv2MatOp(ln,imageSize,F.shape);
% eta = 0.01;
%}
%##### iteration #####
for i = 1 : (iter + 1)
    
    im_residual     =   F * gaussian_dI - im;
    % crop away edge
    kernelSize      =   min(F.xsize, F.fsize);
    corpMarginSize  =   kernelSize;
    Pim             =   patimat('same',size(im_residual),corpMarginSize,0);
    im_residual     =   Pim'*im_residual;
    
    % absolute error
    errorabso       =   gaussian_dI - nature;
    if unique(abs(kernelSize - size(gaussian_dI)) > abs(max(F.xsize, F.fsize) - size(gaussian_dI)))
        errorabso   =   Pim'*errorabso;        % current solving x, crop    
        natureCrop  =   Pim'*nature;
    else
        errorabso   =   errorabso;             % current solving f, NOT crop
        natureCrop  =   nature;
    end
    
    % average residual & relative error of ground truth
    if norm(im_residual,'fro' )==0
        errs(i)     =   1e-20;
        rerrs(i)    =   1e-20;
    else            
        errs(i)     =   (norm(im_residual,'fro') / numel(im_residual)); % average, absolute residual
        rerrs(i)    =   (norm(errorabso,'fro') / norm(natureCrop,'fro')); % relative error ||x - hat(x)|| / ||x||
    end
    % plot
    f2 = figure(2); set(gcf,'visible','off');
    subplot(121) 
    imagesc(clip(gaussian_dI,1,0)); axis image, colormap(gray)
    title(sprintf('Gaussian - Iteration %d/%d',i,iter + 1))
    drawnow    
    subplot(122)
    hData           =   loglog(errs, 'Color', blu);
    hYLabel         =   ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
    hXLabel         =   xlabel('$\#steps$', 'Interpreter','Latex');
    thisFigure;         drawnow 
        
    % stop creterien
    if i > 3
        errRelChange    =  errs(2:i) - errs(1:i-1);
        errRelChange    =  sum(errRelChange(i-3:i-1)) / sum(errs(i-3:i)) * 4/3 ;
    end
    if i == (iter + 1) || errRelChange > inf %-1e-3
        gaussian_dI     =   clip(gaussian_dI,1,0);
        break
    end    
    
    tStart = tic; 
%%%%%%%%%%%%%%%%%%%%%
    % pure
%     gaussian_dI = gaussian_dI .* ( ( F' * im + epsl) ./ ( F'*(F*gaussian_dI) + epsl ) ); % <-- update rule
    
    % Laplacian regularization    
    bgau            =   clip(F'*im, 1e300, 1e-7);
    agau            =   clip( ( F'*(F*gaussian_dI) ) , 1e300, 1e-7) + eta*(lap(gaussian_dI,'+'));
    cgau            =   clip(eta*(lap(gaussian_dI,'-')), inf, -inf);
    
    num             =   clip(bgau + sqrt(clip(bgau.*bgau + 4*agau.*cgau,1e300,1e-24)), 1e300, 1e-7); % + epsl;
    denom           =   2.*agau + epsl;
    update          =   num ./ denom;
    
%     keyboard
    gaussian_dI     =   clip(gaussian_dI .* update, inf, 0);
    
    figure(33), 
    subplot(221), imagesc(update), title('update'), colormap gray, axis image off
    subplot(223), imagesc(gaussian_dI), title('gaussian esti.'), colormap gray, axis image off
    subplot(222), imagesc(num), title('num'), colormap gray, axis image off
    subplot(224), imagesc(denom), title('denom'), colormap gray, axis image off
%%%%%%%%%%%%%%%%%%%%%
    tElapsed = toc(tStart);
    time = [time; time(end)+tElapsed];
end


%##### figure #####
%----- main curves -----
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
if isfield('option','plotFlag') && option.plotFlag == 1
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121), hData = loglog(time ,errs,'Color',blu); thisFigure; hold on
subplot(122), hData = loglog(time,rerrs,'Color',blu); thisFigure; hold on
fstp = figure(15); set(fstp,'visible','on'),
subplot(121), hData = loglog( errs,'Color',blu); thisFigure; hold on
subplot(122), hData = loglog(rerrs,'Color',blu); thisFigure; hold on
% for latex
f10=figure(10); set(f10,'visible','off');
hData = loglog(time, errs,'Color',blu); 
axis tight; thisFigure; hold on
f12=figure(12); set(f12,'visible','off');
hData = loglog(time,rerrs,'Color',blu); 
axis tight; thisFigure; hold on
f11=figure(11); set(f11,'visible','off');
hData = loglog(errs, 'Color',blu); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
f13=figure(13); set(f13,'visible','off');
hData = loglog(rerrs,'Color',blu); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
end
%
%----- image evolution and residual curve -----
figPath = option.figPath;
%
f2 = figure(2); set(f2,'visible','off')
filename = 'deconv_gaussian_with_curve';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% keyboard
%----- gaussian deconved image -----
f_gaussian = figure; set(f_gaussian,'visible','off');
imagesc(clip(gaussian_dI,1,0)); axis image off, colormap(gray)
title('gaussian')
filename = 'deconv_gaussian';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf
end