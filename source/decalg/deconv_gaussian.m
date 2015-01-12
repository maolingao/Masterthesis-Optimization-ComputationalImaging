function [gaussian_dI,errs] = deconv_gaussian(F,im,iter,nature,start,eta,option)
% Gaussian noise deconvolution
startup;

if nargin < 5
    start = F'*im;
    % start = start./sum(start(:)); % nfactor
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
imageSize = size(gaussian_dI);
errs = nan(1,iter);
rerrs = nan(1,iter);
            
epsl = 1e-7;
time = 1e-2;
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
    
    im_residual = F * gaussian_dI - im;
%     im_residual = (gaussian_dI - double(nature)) ;
    if norm(im_residual )==0
        errs(i) = 1e-20;
        rerrs(i) = 1e-20;
    else
        errs(i) = norm(im_residual)/ numel(im_residual); % average, absolute residual
        rerrs(i) = norm(gaussian_dI - nature)/ norm(nature); % relative error ||x - hat(x)|| / ||x||
    end
    
    f2 = figure(2);
    subplot(121) 
    imagesc(clip(gaussian_dI,1,0)); axis image, colormap(gray)
    title(sprintf('Gaussian noise - Iteration %d/%d',i,iter + 1))
    drawnow    
    
    subplot(122)
    loglog(errs,'Color',blu)        
    ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex')
    xlabel('$\#steps$','Interpreter','Latex')
    drawnow
        
    if i == (iter + 1)
        break
    end    
    tStart = tic; 
%%%%%%%%%%%%%%%%%%%%%
    % pure
%     gaussian_dI = gaussian_dI .* ( ( F' * im + epsl) ./ ( F'*(F*gaussian_dI) + epsl ) ); % <-- trivial update rule
    
    % Laplacian regularization    
    bgau = clip(F'*im,inf, 0);
%     agau = clip( ( F'*(F*gaussian_dI) + eta*((Lp*gaussian_dI)) ) ,inf,0);
%     cgau = clip(eta*((Ln*gaussian_dI)), inf, 0);
    agau = clip( ( F'*(F*gaussian_dI) + eta*(lap(gaussian_dI,'+')) ) ,inf,0);
    cgau = clip(eta*(lap(gaussian_dI,'-')), inf, 0);
    
    nom = bgau + sqrt(bgau.*bgau + 4*agau.*cgau) + epsl;
    denom = 2.*agau + epsl;
    
    update = nom ./ denom;
    gaussian_dI = clip(gaussian_dI .* update,1,0);
%%%%%%%%%%%%%%%%%%%%%
    tElapsed = toc(tStart);
    time = [time; time(end)+tElapsed];
end


%##### figure #####
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
%----- main curves -----
% for debug
fclk = figure(14);  set(fclk,'visible','on'),subplot(121),loglog(time,errs,'Color',blu),hold on, 
subplot(122), set(fclk,'visible','on'),loglog(time,rerrs,'Color',blu),hold on
fstp = figure(15); set(fstp,'visible','on'),subplot(121),loglog(errs,'Color',blu),hold on, 
subplot(122), set(fstp,'visible','on'),loglog(rerrs,'Color',blu),hold on
% for latex
f10=figure(10); set(f10,'visible','off'),loglog(time,errs,'Color',blu),hold on, 
f12=figure(12); set(f12,'visible','off'),loglog(time,rerrs,'Color',blu),hold on
f11=figure(11); set(f11,'visible','off'),loglog(errs,'Color',blu),hold on, 
f13=figure(13); set(f13,'visible','off'),loglog(rerrs,'Color',blu),hold on
%
%----- image evolution and residual curve -----
figPath = option.figPath;
%
f2 = figure(2); set(f2,'visible','on')
filename = 'deconv_gaussian_with_curve';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
%----- gaussian deconved image -----
f_gaussian = figure; set(f_gaussian,'visible','off');
imagesc(clip(gaussian_dI,1,0)); axis equal off, colormap(gray)
title('my gaussian')
filename = 'deconv_gaussian';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
end