function [lucy_dI,errs] = deconv_rl(F,im,iter,nature,start,eta,option)
% Richardson-Lucy deconvolution
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

lucy_dI = start; % custermized start guess; 
imageSize = size(lucy_dI);
%##### Tikhonov #####
%{
l = [0 -1 0
     -1 4 -1
     0 -1 0]; % laplacian matrix WORK! 
L = conv2MatOp(l,imageSize,'same');
%}
% eta = 0.001;  % choice according to image noise
% namda =  1e-3;
% ##### setup #####
errs = nan(1,iter);
rerrs = nan(1,iter);
epsl = 1e-7;
time = 1e-2;
for i = 1 : (iter + 1)    

    im_residual = (F * lucy_dI - im) ;
%     im_residual = (lucy_dI - double(nature)) ;
    

    if norm(im_residual ,'fro')==0
        errs(i) = 1e-20;
        rerrs(i) = 1e-20;
    else
        errs(i) = norm(im_residual,'fro')/ numel(im_residual); % average, absolute residual
        rerrs(i) = norm((lucy_dI - nature),'fro')/ norm(nature,'fro'); % relative error ||x - hat(x)|| / ||x||
    end

    
    f1 = figure(1);
    subplot(121)
    imagesc(clip(lucy_dI,1,0)); axis image, colormap(gray)
    title(sprintf('lucy - Iteration %d/%d',i,iter + 1))
    drawnow    
    
    subplot(122)
    hData = loglog(errs, 'Color', ora);
    hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
    hXLabel = xlabel('$\#steps$', 'Interpreter','Latex');
    thisFigure;   
    drawnow 
    
    if i == (iter + 1)
        break
    end
    if i == 50
        NOP
    end

    tStart = tic; 
%%%%%%%%%%%%%%%%%%%%%
    % ***** poisson noise model *****
    % pure
%     lucy_dI = lucy_dI .* ( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl);
    % TV regularization   
    %{
    [gx,gy] = gradient(lucy_dI);
    nablaI = sqrt(gx.^2 + gy.^2);
    [dgx, ~] = gradient(gx./(nablaI + epsl));
    [~, dgy] = gradient(gy./(nablaI + epsl));
    rt = - (dgx + dgy); 
    lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* ( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl);
    lucy_dI = clip(lucy_dI,1,0);
    %}
    % Laplacian regularization
    %{
    LI = clip(laplacian(lucy_dI), inf, -inf);
    ItrLI = clip(lucy_dI'*LI, inf, -inf);                          % Itranspose*laplacian*I
    [gx,gy] = gradient(lucy_dI);
    nablaI = sqrt(gx.^2 + gy.^2);
    etaf = -1/eta;
%     rt = clip(etaf * exp( etaf * ItrLI) * nablaI * LI, inf, -inf); % Itranspose*laplacian*I
%     rt = clip(etaf * exp( etaf * nablaI.*nablaI) * nablaI * LI, inf, -inf); % abs(nabla*I)
%     rt = clip( exp( ItrLI) .* nablaI .* LI, inf, -inf);               % Itranspose*laplacian*I
    rt = clip( exp( nablaI.*nablaI) .* nablaI .* LI, inf, -inf); % abs(nabla*I)
%     lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* ( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl);
    lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* (( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl));
    lucy_dI = clip(lucy_dI,1,0);
    %}
    %{
    % Tai Paper
    LI = lap(lucy_dI);                          % laplace(I)
    [gx,gy] = gradient(lucy_dI);
    nablaI = sqrt(gx.^2 + gy.^2);                     % abs(nabla(I))
    rt = exp( nablaI.*nablaI) .* nablaI .* LI;
    lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* (( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl));         % clip pixel here btw 0 and 1
    lucy_dI = clip(lucy_dI,1,0);    
    %}
    % Tai code    
    [gx,gy] = gradient(lucy_dI);
    [ggx,~] = gradient(flip(gx,2));
    [~,ggy] = gradient(flip(gy,2));
    rt = abs(gx).*ggx + abs(gy).*ggy ;
    
     lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* (( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl));         % clip pixel here btw 0 and 1
%    lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* ( F' * ((im ) ./ (F*lucy_dI + epsl) )  );         % clip pixel here btw 0 and 1
   
    lucy_dI = clip(lucy_dI,1,0);    
    %{
    % easy update by Levin
    LI = clip(laplacian(lucy_dI), inf, -inf);
    rt = LI;
    lucy_dI = lucy_dI ./ (1 - eta * rt + epsl) .* (( F' * ((im  + epsl ) ./ (F*lucy_dI + epsl) ) + epsl  ) ./ ( F' * ones(size(im)) + epsl));
    lucy_dI = clip(lucy_dI,1,0);
    %}
    % ***** gaussian noise model *****
    % pure
%     lucy_dI = lucy_dI + ( F' * (im  - (F * lucy_dI) ) );
    % TV regularization    
    %{
    [gx,gy] = gradient(lucy_dI);
    nablaI = sqrt(gx.^2 + gy.^2);
    [dgx, ~] = gradient(gx./(nablaI + epsl));
    [~, dgy] = gradient(gy./(nablaI + epsl));
    rt = - (dgx + dgy); 
    lucy_dI = lucy_dI + ( F' * (im  - (F * lucy_dI) ) ) + eta * rt;
    lucy_dI = clip(lucy_dI,1,0);
    %}
    % Laplacian regularization
    %{
    LI = clip(laplacian(lucy_dI), inf, -inf);
    ItrLI = clip(lucy_dI'*LI, inf, -inf);                          % Itranspose*laplacian*I
    [gx,gy] = gradient(lucy_dI);
    nablaI = sqrt(gx.^2 + gy.^2);
    etaf = -1/eta;
%     rt = clip(etaf * exp( etaf * ItrLI) .* nablaI .* LI, 1, -1); % Itranspose*laplacian*I
%     rt = clip(etaf * exp( etaf * nablaI.*nablaI) * nablaI * LI, inf, -inf); % abs(nabla*I)
    rt = clip( exp( ItrLI) .* nablaI .* LI, 1, -1); % Itranspose*laplacian*I
%     rt = clip( exp(  nablaI.*nablaI) .* nablaI .* LI, inf, -inf); % abs(nabla*I)
    lucy_dI = lucy_dI + ( F' * (im  - (F * lucy_dI) ) ) + eta * rt;
    lucy_dI = clip(lucy_dI,1,0);
%     figure(10), imagesc(rt), colormap gray, axis equal
    %}
    % Bilateral regularization
%%%%%%%%%%%%%%%%%%%%%
    tElapsed = toc(tStart);
    time = [time;time(end)+tElapsed];

end

%%##### figure #####
%----- main curves -----
errs = errs(~isnan(errs));
rerrs = rerrs(~isnan(rerrs));
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121), hData = loglog(time ,errs,'Color',ora); thisFigure; hold on
subplot(122), hData = loglog(time,rerrs,'Color',ora); thisFigure; hold on
fstp = figure(15); set(fstp,'visible','on'),
subplot(121), hData = loglog( errs,'Color',ora); thisFigure; hold on
subplot(122), hData = loglog(rerrs,'Color',ora); thisFigure; hold on
% for latex
f10=figure(10); set(f10,'visible','off');
hData = loglog(time, errs,'Color',ora); 
axis tight; thisFigure; hold on
f12=figure(12); set(f12,'visible','off');
hData = loglog(time,rerrs,'Color',ora); 
axis tight; thisFigure; hold on
f11=figure(11); set(f11,'visible','off');
hData = plot(errs, 'Color',ora); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
f13=figure(13); set(f13,'visible','off');
hData = plot(rerrs,'Color',ora); 
set(gca,'Yscale','log'), axis tight; thisFigure; hold on 
%
%----- image evolution and residual curve -----
figPath = option.figPath;

f1 = figure(1); set(f1,'visible','on')
filename = 'deconv_lucy_with_curve';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% keyboard
%----- lucy deconved image -----
f_lucy = figure; set(f_lucy,'visible','off');
imagesc(clip(lucy_dI,1,0)); axis image off, colormap(gray)
title('lucy')
filename = 'deconv_lucy';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
end