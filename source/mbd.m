function I = mbd(multiFrame, F, start, iterK, natureI, multiKernel, option)
% [pncg_dI, cg_dI, lucy_dI, gaussian_dI] = mbd(multiFrame, F, start, iterK, natureI, multiKernel)
% offline multi-frame blind deconvolution
if nargin < 7
    option.figPath = '/is/ei/mgao/figure2drag';
    option.method = 'gaussian';
end

startup;
figPath = option.figPath;

f10 = figure(10); clf(10), set(f10,'visible','off')
f11 = figure(11); clf(11), set(f11,'visible','off')
f12 = figure(12); clf(12), set(f12,'visible','off')
f13 = figure(13); clf(13), set(f13,'visible','off')
fclk = figure(14); clf(fclk), set(fclk,'visible','on')
fstp = figure(15); clf(fstp), set(fstp,'visible','on')

if F.xsize > F.fsize
    imagesize = F.xsize;
    fsize = F.fsize;
else
    imagesize = F.fsize;
    fsize = F.xsize;
end
shape = F.shape;

errs_allframes_pncg = []; rerrs_allframes_pncg = [];
errs_allframes_cg = []; rerrs_allframes_cg = [];
errs_allframes_gaussian = []; rerrs_allframes_gaussian = [];
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
eta = 0;                                        % ### <--- regularization parameter
% kernel estimating
tolK = -inf;
scaler = 1e0;
HK = hessianMatrix(eye(fsize)*scaler);
X = conv2MatOp(im2double(start),fsize,shape);   % initial guess of convMtx X
% -------- ground truth comparison figure - start --------
startImg = figure;  set(startImg,'visible','off'),
subplot(1,2,1)
imagesc(clip(natureI,1,0)); 
axis image off; colormap gray;
subplot(1,2,2)
imagesc(clip(start,1,0)); 
axis image off; colormap gray;
filename = sprintf('startNatureImg_%d',0);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
% ############### gradient img to deconvolve f ###############
% keyboard
% gstart = lap(start);
% % gstart = gstart - min(vec(gstart));
% % gstart = gstart./max(vec(gstart));
% % figure, imagesc(gstart); colormap gray, axis image
% X = conv2MatOp(im2double(gstart),fsize,shape);   % initial guess of convMtx X
% ############### END of gradient img to deconvolve f ###############


% nature estimating
iterN = 1; % one step for estimating nature
HN = hessianMatrix(eye(imagesize));
tolN = -inf;
% iteration
%{
% all methods in one loop
for k = 1 : numFrame
    % ##### special setup #####
    frame = multiFrame{k};
    natureK = multiKernel{k}; % for error calculation
    option.version = 'FH';
    % ##### estimate kernel #####
    [pncg_kernel,HK,errs_pncg_kernel] = deconv_pncg(X,frame,natureK,HK,iterK,startK,tolK,eta,option); % pncg
    [cg_kernel,errs_cg_kernel] = deconv_cg(X, frame, natureK, iterK, startK, tolK, eta, option); % cg
    [gaussian_kernel,errs_gaussian_kernel] = deconv_gaussian(X,frame,iterK,natureK,startK,eta); % gaussian
    % ##### estimate nature #####
    Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
    Kcg = conv2MatOp(im2double(cg_kernel),imagesize,shape);  % convMtx of kernel, cg
    Kgaussian = conv2MatOp(im2double(gaussian_kernel),imagesize,shape);  % convMtx of kernel, gaussian
    if k == 1
        [pncg_dI,~,errs_pncgN] = deconv_pncg(Kpncg,frame,natureI,HN,iterN,start,tolN,eta,option); % pncg
        [cg_dI, errs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, start, tolN, eta, option); % cg
        [gaussian_dI,errs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,start,eta); % gaussian
    else
        [pncg_dI,~,errs_pncgN] = deconv_pncg(Kpncg,frame,natureI,HN,iterN,pncg_dI,tolN,eta,option); % pncg
        [cg_dI,errs_cgN] = deconv_cg(Kcg, frame, natureI, HN, iterN, cg_dI, tolN, eta, option); % cg
        [gaussian_dI,errs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,gaussian_dI,eta); % gaussian
    end
    
    % statitics 
    % all frame errors
    errs_allframes_pncg = [errs_allframes_pncg,errs_pncgN(end)]; % pncg 
    errs_allframes_cg = [errs_allframes_cg,errs_cgN(end)]; % cg 
    errs_allframes_gaussian = [errs_allframes_gaussian,errs_gaussianN(end)]; % gaussian 
    
    % plots    
    % pncg
    figure(f_pncg); subplot(121)
    imagesc(clip(pncg_dI,1,0)), axis image,colormap(gray)
    title(sprintf('pncg - frame %d/%d',k,length(multiFrame)))
    drawnow
    subplot(122)   
    plot(1:length(errs_allframes_pncg),errs_allframes_pncg,'k')
    xlabel('#frames'), ylabel('residual |Fx - y| / pixel')
    title(sprintf('frame %d/%d',k,numFrame))
    drawnow    
    % cg
    figure(f_cg); subplot(121)
    imagesc(clip(cg_dI,1,0)), axis image,colormap(gray)
    title(sprintf('cg - frame %d/%d',k,length(multiFrame)))
    drawnow
    subplot(122)   
    plot(1:length(errs_allframes_cg),errs_allframes_cg,'k')
    xlabel('#frames'), ylabel('residual |Fx - y| / pixel')
    title(sprintf('frame %d/%d',k,numFrame))
    drawnow    
    % gaussian
    figure(f_gaussian); subplot(121)
    imagesc(clip(gaussian_dI,1,0)), axis image,colormap(gray)
    title(sprintf('gaussian - frame %d/%d',k,length(multiFrame)))
    drawnow
    subplot(122)   
    plot(1:length(errs_allframes_gaussian),errs_allframes_gaussian,'k')
    xlabel('#frames'), ylabel('residual |Fx - y| / pixel')
    title(sprintf('frame %d/%d',k,numFrame))
    drawnow    
end
%}
switch option.method
    case 'pncg'
        startK = zeros(fsize);                   % initial guess of kernel, flat image with all NULL
        f_pncg = figure(114); clf(f_pncg); set(f_pncg,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            option.version = 'FH';
            % ##### estimate kernel #####
            % !!!!!!!! scale down image deconvolution matrix X and blurry observation frame !!!!!!!!
%             if k == 1
%                 natureI = natureI./10^(3);   % <---- scale X
%             end
%             frame = frame./10^(3);           % <---- scale y
            % !!!!!!!! non blind !!!!!!!!
            clear X
            X = conv2MatOp(im2double(natureI),fsize,shape);   % convMtx X of nature --- mfd
%             clear HK
%             HK = hessianMatrix(eye(fsize)*1e0);               % same initial for HessianMatrix HK --- mfd
            %
            [pncg_kernel, HK, errs_pncg_kernel] = deconv_pncg(X, frame, natureK, HK, iterK, startK, tolK, eta, option); % pncg
            pncg_kernel = preserveNorm(pncg_kernel);            % preserve energy norm of PSF
            
            % ##### MEMLIM #####
            %{%
            MEMLIM = 50;% size(HK.s,2);
%             keyboard
            [S,Y,Delta,GInv] = purify(HK.s,HK.y,HK.delta,HK.Ginv0,MEMLIM);
            % keyboard
            clear HK
            HK = hessianMatrix(eye(fsize)*scaler, S, Y, Delta, GInv, size(S,2)+1);
%             keyboard
            max(max(S'*Y - diag(1./diag(HK.Ginv0))))
            %}
            % ###################
            % -------- kernel comparison figure --------
            pncgKernelImg = figure; set(pncgKernelImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureK,1,0)); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(clip(pncg_kernel,1,0)); 
            axis image off; colormap gray;
            filename = sprintf('pncgKernelImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
%             keyboard
            % ##### estimate nature #####
            clear Kpncg
            Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
            if k == 1
                [pncg_dI, ~, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, start, tolN, eta, option); % pncg
            else
                [pncg_dI, ~, errs_pncgN, ~, rerrs_pncgN] = deconv_pncg(Kpncg, frame, natureI, HN, iterN, pncg_dI, tolN, eta, option); % pncg
            end
%             if k == 1
%                 [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,start,eta); % gaussian
%             else
%                 [pncg_dI,errs_pncgN,rerrs_pncgN] = deconv_gaussian(Kpncg,frame,iterN,natureI,pncg_dI,eta); % gaussian
%             end
            pncg_dI = clip(pncg_dI,inf,0);
            clear X
            X = conv2MatOp(im2double(pncg_dI),fsize,shape);   % new guess of convMtx X %#################
            % -------- ground truth comparison figure --------
            pncgNatureImg = figure;  set(pncgNatureImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureI,1,0)); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(clip(pncg_dI,1,0)); 
            axis image off; colormap gray;
            filename = sprintf('pncgNatureImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            % statitics 
            % all frame errors
            errs_allframes_pncg = [errs_allframes_pncg,errs_pncgN(end)]; % pncg 
            rerrs_allframes_pncg = [rerrs_allframes_pncg,rerrs_pncgN(end)]; % relative error, pncg

            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            figure(f_pncg); subplot(121)
            imagesc(clip(pncg_dI,1,0)), axis image,colormap(gray)
            title(sprintf('pncg - frame %d/%d',k,length(multiFrame)))
            drawnow
            subplot(122)   
            hData = plot(1:length(errs_allframes_pncg),errs_allframes_pncg,'Color', dre);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('frame %d/%d',k,numFrame));
            thisFigure;   
            drawnow        
        end
            % -------- ground truth frame error figure --------
            % for latex
            % residual error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(errs_allframes_pncg),errs_allframes_pncg,'Color', dre);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_errAllFrame_pncg');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            % relative error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(rerrs_allframes_pncg),rerrs_allframes_pncg,'Color', dre);          
            hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_relativeErrorAllFrame_pncg');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
        I = pncg_dI;
    case 'cg'
        startK = zeros(fsize);                   % initial guess of kernel, flat image with all NULL
        f_cg = figure(115); clf(f_cg); set(f_cg,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            % ############### gradient img to deconvolve f ###############
%             keyboard
%             gframe = lap(frame);
%             gframe = gframe - min(vec(gframe));
%             gframe = gframe./max(vec(gframe));
%             figure, imagesc(gframe); colormap gray, axis image
            % ############### END of gradient img to deconvolve f ###############
            %
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            % !!!!!!!! non blind !!!!!!!!
            clear X
            X = conv2MatOp(im2double(natureI),fsize,shape);   % convMtx X of nature --- mfd
            %
            %
            [cg_kernel,errs_cg_kernel] = deconv_cg(X, frame, natureK, iterK, startK, tolK, eta, option); % cg
            cg_kernel = preserveNorm(cg_kernel);            % preserve energy norm of PSF
            % -------- kernel comparison figure --------
            cgKernelImg = figure; set(cgKernelImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureK,1,0)); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(clip(cg_kernel,1,0)); 
            axis image off; colormap gray;
            filename = sprintf('cgKernelImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
%             keyboard
            % ##### estimate nature #####
            clear Kcg
            Kcg = conv2MatOp(im2double(cg_kernel),imagesize,shape);  % convMtx of kernel, cg
            if k == 1
                [cg_dI, errs_cgN, ~, rerrs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, start, tolN, eta, option); % cg
            else
                [cg_dI, errs_cgN, ~, rerrs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, cg_dI, tolN, eta, option); % cg
            end
%             if k == 1
%                 [cg_dI,errs_cgN] = deconv_gaussian(Kcg,frame,iterN,natureI,start,eta); % gaussian
%             else
%                 [cg_dI,errs_cgN] = deconv_gaussian(Kcg,frame,iterN,natureI,cg_dI,eta); % gaussian
%             end
%             keyboard
%             cg_dI = feasible(cg_dI);                        % scale and clip pixel value [0,1]
            cg_dI = clip(cg_dI,1,0);
            clear X
            X = conv2MatOp(im2double(cg_dI),fsize,shape);   % new guess of convMtx X
            % ############### gradient img to deconvolve f ###############
%             gcg_dI = lap(cg_dI);
% %             gcg_dI = gcg_dI - min(vec(gcg_dI));
% %             gcg_dI = gcg_dI./max(vec(gcg_dI));
%             X = conv2MatOp(im2double(gcg_dI),fsize,shape);   % new guess of convMtx X
            % ############### END of gradient img to deconvolve f ###############
            %            
            %--------------- TOOL ---------------%
            % sum of cg_dI(more an edge img) and frame
            %
%             keyboard
%             cg_dI_show = bsxfun(@times,cg_dI,(cg_dI>1e-5));
%             cg_dI_show = bsxfun(@times, cg_dI_show, sum(vec(frame))/sum(vec(cg_dI_show)));
%             cg_dI_show = cg_dI_show + frame;
%             cg_dI_show = bsxfun(@times, cg_dI_show, sum(vec(frame))/sum(vec(cg_dI_show)));
            %--------------- END  ---------------%
            % -------- ground truth comparison figure --------
            cgNatureImg = figure;  set(cgNatureImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureI,1,0)); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(clip(cg_dI,1,0)); 
            axis image off; colormap gray;
            filename = sprintf('cgNatureImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;

            % statitics 
            % all frame errors 
            errs_allframes_cg = [errs_allframes_cg,errs_cgN(end)]; % residual, cg 
            rerrs_allframes_cg = [rerrs_allframes_cg,rerrs_cgN(end)]; % relative error, cg
            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            figure(f_cg); subplot(121)
            imagesc(clip(cg_dI,1,0)), axis image,colormap(gray)
            title(sprintf('cg - frame %d/%d',k,length(multiFrame)))
            drawnow
            subplot(122)   
            hData = plot(1:length(errs_allframes_cg),errs_allframes_cg,'Color', mpg);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('frame %d/%d',k,numFrame));
            thisFigure;   
            drawnow        
        end 
            % -------- ground truth frame error figure --------
            % for latex
            % residual error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(errs_allframes_cg),errs_allframes_cg,'Color', mpg);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_residualErrorAllFrame_cg');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            % relative error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(rerrs_allframes_cg),rerrs_allframes_cg,'Color', mpg);          
            hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_relativeErrorAllFrame_cg');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
        I = cg_dI;
    case 'rl'
    case 'gaussian'
        startK = ones(fsize)./prod(fsize);              % initial guess of kernel, flat image with uniform entries
        f_gaussian = figure(112); clf(f_gaussian); set(f_gaussian,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            % !!!!!!!! non blind !!!!!!!!
            clear X
            X = conv2MatOp(im2double(natureI),fsize,shape);   % convMtx X of nature --- mfd
            %
            [gaussian_kernel,errs_gaussian_kernel] = deconv_gaussian(X,frame,iterK,natureK,startK,eta,option); % gaussian
            gaussian_kernel = preserveNorm(gaussian_kernel);            % preserve energy norm of PSF
            % -------- kernel comparison figure --------
            gaussianKernelImg = figure; set(gaussianKernelImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureK,1,0)); 
            axis image off; colormap gray;
            drawnow
            subplot(1,2,2)
            imagesc(clip(gaussian_kernel,1,0)); 
            axis image off; colormap gray;
            drawnow
            filename = sprintf('gaussianKernelImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
%             keyboard
            % ##### estimate nature #####
            clear Kgaussian
            Kgaussian = conv2MatOp(im2double(gaussian_kernel),imagesize,shape);  % convMtx of kernel, gaussian
            if k == 1
                [gaussian_dI,errs_gaussianN,rerrs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,start,eta); % gaussian
            else
                [gaussian_dI,errs_gaussianN,rerrs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,gaussian_dI,eta); % gaussian
            end
            clear X
            X = conv2MatOp(im2double(gaussian_dI),fsize,shape);   % new guess of convMtx X
            % -------- ground truth comparison figure --------
            gaussianNatureImg = figure;  set(gaussianNatureImg,'visible','off'),
            subplot(1,2,1)
            imagesc(clip(natureI,1,0)); 
            axis image off; colormap gray;
            drawnow
            subplot(1,2,2)
            imagesc(clip(gaussian_dI,1,0)); 
            axis image off; colormap gray;
            drawnow
            filename = sprintf('gaussianNatureImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
            % statitics 
            % all frame errors
            errs_allframes_gaussian = [errs_allframes_gaussian,errs_gaussianN(end)]; % gaussian 
            rerrs_allframes_gaussian = [rerrs_allframes_gaussian,rerrs_gaussianN(end)]; % relative error, gaussian

            % plots    
            % -------- ground truth frame error figure --------
            % for debug
            figure(f_gaussian); subplot(121)
            imagesc(clip(gaussian_dI,1,0)), axis image,colormap(gray)
            title(sprintf('gaussian - frame %d/%d',k,length(multiFrame)))
            drawnow
            subplot(122)   
            hData = plot(1:length(errs_allframes_gaussian),errs_allframes_gaussian,'Color', blu);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('frame %d/%d',k,numFrame));
            thisFigure;   
            drawnow
            
        end            
            % -------- ground truth frame error figure --------
            % for latex
            % residual error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(errs_allframes_gaussian),errs_allframes_gaussian,'Color', blu);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_residualErrorAllFrame_gaussian');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            % relative error
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(rerrs_allframes_gaussian),rerrs_allframes_gaussian,'Color', mpg);          
            hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            drawnow
            filename = sprintf('mbd_relativeErrorAllFrame_gaussian');
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
        I = gaussian_dI;
    otherwise
        print '[mbd.m]: option.method can be one of those: 'pncg', 'cg', 'rl', 'gaussian''
end
    

%##### figure #####
saveResultFigure;

end
