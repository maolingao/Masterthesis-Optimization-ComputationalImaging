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

errs_allframes_pncg = [];
errs_allframes_cg = [];
errs_allframes_gaussian = [];
% ##### general setup #####
numFrame = numel(multiFrame);
eta = 0.05;
% kernel estimating
tolK = 1e-7;
HK = hessianMatrix(eye(fsize));
X = conv2MatOp(im2double(start),fsize,shape);   % initial guess of convMtx X
startK = ones(fsize)./prod(fsize);              % initial guess of kernel, flat image
% nature estimating
iterN = 1; % one step for estimating nature
HN = hessianMatrix(eye(imagesize));
tolN = 1e-5;
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
        f_pncg = figure(114); clf(f_pncg); set(f_pncg,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            option.version = 'FH';
            % ##### estimate kernel #####
            [pncg_kernel,HK,errs_pncg_kernel] = deconv_pncg(X,frame,natureK,HK,iterK,startK,tolK,eta,option); % pncg
            % ##### estimate nature #####
            Kpncg = conv2MatOp(im2double(pncg_kernel),imagesize,shape);  % convMtx of kernel, pncg
            if k == 1
                [pncg_dI,~,errs_pncgN] = deconv_pncg(Kpncg,frame,natureI,HN,iterN,start,tolN,eta,option); % pncg
            else
                [pncg_dI,~,errs_pncgN] = deconv_pncg(Kpncg,frame,natureI,HN,iterN,pncg_dI,tolN,eta,option); % pncg
            end

            % statitics 
            % all frame errors
            errs_allframes_pncg = [errs_allframes_pncg,errs_pncgN(end)]; % pncg 

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
        end
        I = pncg_dI;
    case 'cg'
        f_cg = figure(115); clf(f_cg); set(f_cg,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            [cg_kernel,errs_cg_kernel] = deconv_cg(X, frame, natureK, iterK, startK, tolK, eta, option); % cg
            % ##### estimate nature #####
            Kcg = conv2MatOp(im2double(cg_kernel),imagesize,shape);  % convMtx of kernel, cg
            if k == 1
                [cg_dI, errs_cgN] = deconv_cg(Kcg, frame, natureI, iterN, start, tolN, eta, option); % cg
            else
                [cg_dI,errs_cgN] = deconv_cg(Kcg, frame, natureI, HN, iterN, cg_dI, tolN, eta, option); % cg
            end

            % statitics 
            % all frame errors 
            errs_allframes_cg = [errs_allframes_cg,errs_cgN(end)]; % cg 

            % plots    
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
        end
        I = cg_dI;
    case 'rl'
    case 'gaussian'
        f_gaussian = figure(112); clf(f_gaussian); set(f_gaussian,'visible','on')
        for k = 1 : numFrame
            % ##### special setup #####
            frame = multiFrame{k};
            natureK = multiKernel{k}; % for error calculation
            % ##### estimate kernel #####
            [gaussian_kernel,errs_gaussian_kernel] = deconv_gaussian(X,frame,iterK,natureK,startK,eta,option); % gaussian
            % -------- kernel comparison figure --------
            gaussianKernelImg = figure; set(gaussianKernelImg,'visible','off'),
            subplot(1,2,1)
            imagesc(natureK); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(gaussian_kernel); 
            axis image off; colormap gray;
            filename = sprintf('gaussianKernelImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
%             keyboard
            % ##### estimate nature #####
            clear Kgaussian
            Kgaussian = conv2MatOp(im2double(gaussian_kernel),imagesize,shape);  % convMtx of kernel, gaussian
            if k == 1
                [gaussian_dI,errs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,start,eta); % gaussian
            else
                [gaussian_dI,errs_gaussianN] = deconv_gaussian(Kgaussian,frame,iterN,natureI,gaussian_dI,eta); % gaussian
            end
            clear X
            X = conv2MatOp(im2double(gaussian_dI),fsize,shape);   % new guess of convMtx X
            % -------- ground truth comparison figure --------
            gaussianNatureImg = figure;  set(gaussianNatureImg,'visible','off'),
            subplot(1,2,1)
            imagesc(natureI); 
            axis image off; colormap gray;
            subplot(1,2,2)
            imagesc(gaussian_dI); 
            axis image off; colormap gray;
            filename = sprintf('gaussianNatureImg_%d',k);
            filename = fullfile(figPath,filename);
            print(gcf, '-depsc2', filename)
            close gcf;
            
            % statitics 
            % all frame errors
            errs_allframes_gaussian = [errs_allframes_gaussian,errs_gaussianN(end)]; % gaussian 

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
            figure;  set(gcf,'visible','off'),
            hData = plot(1:length(errs_allframes_gaussian),errs_allframes_gaussian,'Color', blu);          
            hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
            hTitle = title(sprintf('%d frames', numFrame));
            thisFigure;   
            filename = sprintf('mbd_errAllFrame_gaussian');
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
