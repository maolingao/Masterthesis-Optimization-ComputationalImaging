function I = mfd(multiFrame,multiFilt,F,iter,nature,option)
% offline multi-frame non-blind deconvolution
if nargin < 6
    option.figPath = '/is/ei/mgao/figure2drag';
    option.method = 'gaussian';
    option.version = 'FH';
end
figPath = option.figPath;
startup;
f10 = figure(10); clf(10), set(f10,'visible','off')
f11 = figure(11); clf(11), set(f11,'visible','off')
f12 = figure(12); clf(12), set(f12,'visible','off')
f13 = figure(13); clf(13), set(f13,'visible','off')
fclk = figure(14); clf(fclk), set(fclk,'visible','on')
fstp = figure(15); clf(fstp), set(fstp,'visible','on')

if F.xsize>F.fsize
    imagesize = F.xsize;
else
    imagesize = F.fsize;
end
% nature = 0;
errs_allsteps_pncg = [];
errs_allsteps_cg = [];
errs_allsteps_lucy = [];
errs_allsteps_gauss = [];
errs_allframes_pncg = [];
errs_allframes_cg = [];
errs_allframes_lucy = [];
errs_allframes_gauss = [];

% ##### general setup #####
shape = F.shape;
f = multiFilt{1};
fsize = size(f);
switch option.target
    case 'kernel'
        F = conv2MatOp(im2double(nature),fsize,shape);  % <-- guess kernel
    case 'nature'
        NOP
    otherwise
        print '[mfd.m]: option.target can be one of those: 'kernel', 'nature''
end
H = hessianMatrix(eye(fsize));
tol = 1*10^(-30);
eta = 0;
% iteration
for k = 1 : numel(option.method)
    switch option.method{k}
        case 'pncg'
            for i = 1 : length(multiFrame)
                % ##### special setup #####
                switch option.target
                    case 'kernel'
                        nature = multiFilt{i};                          % <-- guess kernel
                        frame = multiFrame{i};
                        start = F'*frame; start = start./sum(start(:)); % nfactor
                        tgtSize = fsize;
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    case 'nature'
                        f = multiFilt{i};                               % <-- guess ground truth
                        F = conv2MatOp(f,imagesize,shape);              % <-- guess ground truth
                        frame = multiFrame{i};
                        start = F'*frame; start = start./max(start(:)); % nfactor
                        tgtSize = size(start);
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    otherwise
                        print '[mfd.m]: option.target can be one of those: 'kernel', 'nature''
                end
                % ##### estimate nature #####
                [pncg_dI,H,data_pncg] = deconv_pncg(F,frame,nature,H,iter,start,tol,eta,option); % pncg
                % ########### MEMLIM #############
                MEMLIM = option.MEMLIM;
                lambda = option.MEMSTR;
                % ------------------------------ %
                % ### evd
                [S,Y,Delta,GInv] = purify(H.s,H.y,H.delta,MEMLIM,lambda);
                clear H
                H = hessianMatrix(eye(tgtSize), S, Y, Delta, [],[]);
                % ------------------------------ %
                % statitics        
                % all step errors
                errs_allsteps_pncg =  [errs_allsteps_pncg, data_pncg.errs]; % pncg
                % all frame errors
                errs_allframes_pncg = [errs_allframes_pncg,data_pncg.errs(end)]; % pncg 
            %
            % plots              
                % residual curves
                figure(200), 
                plot(0:length(data_pncg.errs)-1,data_pncg.errs); hold on
                %----- pncg deconved image -----
                f_pncg = figure; set(f_pncg,'visible','off');
                imagesc(clip(pncg_dI,1,0)); axis image, colormap(gray)
                title('pncg')
                filename = sprintf('deconv_pncg_%d_frame',i);
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)

            end
            I = pncg_dI;        
        case 'cg' 
            for i = 1 : length(multiFrame)
                % ##### special setup #####
                switch option.target
                    case 'kernel'
                        nature = multiFilt{i};                          % <-- guess kernel
                        frame = multiFrame{i};
                        start = F'*frame; start = start./sum(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    case 'nature'
                        f = multiFilt{i};                               % <-- guess ground truth
                        F = conv2MatOp(f,imagesize,shape);              % <-- guess ground truth
                        frame = multiFrame{i};
                        start = F'*frame; start = start./max(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    otherwise
                        print '[mfd.m]: option.target can be one of those: 'kernel', 'nature''
                end
%                 start = zeros(size(F'*frame));
                % ##### estimate nature #####
                if i == 1       
                    [cg_dI,errs_cg] = deconv_cg(F,frame,nature,iter,start,tol,eta,option); % cg
                else
                    [cg_dI,errs_cg] = deconv_cg(F,frame,nature,iter,start,tol,eta,option); % cg
                end
            %
            % statitics        
                % all step errors
                errs_allsteps_cg =  [errs_allsteps_cg, errs_cg]; % cg
                % all frame errors
                errs_allframes_cg = [errs_allframes_cg,errs_cg(end)]; % cg 
            %
            % plots    
                % cg
                f_cg = figure; set(f_cg,'visible','off')
                subplot(121)
                imagesc(clip(cg_dI,1,0)); axis image,colormap(gray)
                title(sprintf('cg - frame %d/%d',i,length(multiFrame)))
                drawnow       
                subplot(122)
                plot(1:length(errs_allframes_cg),errs_allframes_cg,'Color',mpg)
                xlabel('#frames','FontSize', 15), ylabel('residual |Fx - y| / pixel','FontSize', 15)
                title(sprintf('frame %d/%d',i,length(multiFrame)))
                drawnow    
                % save figure
                figure(f_cg); set(f_cg,'visible','off')
                filename = ['mfd_deconv_cg_with_curve'];
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)
                %----- cg deconved image -----
                f_cg = figure; set(f_cg,'visible','off');
                imagesc(clip(cg_dI,1,0)); axis image, colormap(gray)
                title('cg')
                filename = sprintf('deconv_cg_%d_frame',i);
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)

            end
            I = cg_dI;
        case 'rl' 
            for i = 1 : length(multiFrame)
                % ##### special setup #####
                switch option.target
                    case 'kernel'
                        nature = multiFilt{i};                          % <-- guess kernel
                        frame = multiFrame{i};
                        start = F'*frame; start = start./sum(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    case 'nature'
                        f = multiFilt{i};                               % <-- guess ground truth
                        F = conv2MatOp(f,imagesize,shape);              % <-- guess ground truth
                        frame = multiFrame{i};
                        start = F'*frame; start = start./max(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    otherwise
                        print '[mfd.m]: option.target can be one of those: 'kernel', 'nature''
                end
                % ##### estimate nature #####
                if i == 1       
                    [lucy_dI,errs_lucy] = deconv_rl(F,frame,iter,nature,start,eta,option); % lucy
                else
                    [lucy_dI,errs_lucy] = deconv_rl(F,frame,iter,nature,start,eta,option); % lucy
                end
            %
            % statitics        
                % all step errors
                errs_allsteps_lucy =  [errs_allsteps_lucy, errs_lucy]; % lucy
                % all frame errors
                errs_allframes_lucy = [errs_allframes_lucy,errs_lucy(end)]; % lucy
            %
            % plots    
                % lucy
                f_rl = figure;  set(f_rl,'visible','off')
                subplot(121)
                imagesc(clip(lucy_dI,1,0)); axis image, colormap(gray)
                title(sprintf('lucy - frame %d/%d',i,length(multiFrame)))
                drawnow                  
                subplot(122)
                plot(1:length(errs_allframes_lucy),errs_allframes_lucy,'Color',ora)
                xlabel('#frames','FontSize', 15), ylabel('residual |Fx - y| / pixel','FontSize', 15)
                title(sprintf('frame %d/%d',i,length(multiFrame)))
                drawnow    
                % save figure
                figure(f_rl); set(f_rl,'visible','off')
                filename = ['mfd_deconv_rl_with_curve'];
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)
                %----- rl deconved image -----
                f_pncg = figure; set(f_pncg,'visible','off');
                imagesc(clip(lucy_dI,1,0)); axis image, colormap(gray)
                title('lucy')
                filename = sprintf('deconv_lucy_%d_frame',i);
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)

            end
            I = lucy_dI;
        case 'gaussian' 
            for i = 1 : length(multiFrame)
                % ##### special setup #####
                switch option.target
                    case 'kernel'
                        nature = multiFilt{i};                          % <-- guess kernel
                        frame = multiFrame{i};
                        start = F'*frame; start = start./sum(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    case 'nature'
                        f = multiFilt{i};                               % <-- guess ground truth
                        F = conv2MatOp(f,imagesize,shape);              % <-- guess ground truth
                        frame = multiFrame{i};
                        start = F'*frame; start = start./max(start(:)); % nfactor
                        %     start = ones(size(F'*frame))/numel(F'*frame);
                    otherwise
                        print '[mfd.m]: option.target can be one of those: 'kernel', 'nature''
                end
                % ##### estimate nature #####
                if i == 1       
                    [gaussian_dI,errs_gauss] = deconv_gaussian(F,frame,iter,nature,start,eta,option); % gauss
                else
                    [gaussian_dI,errs_gauss] = deconv_gaussian(F,frame,iter,nature,start,eta,option); % gauss
                end
            %
            % statitics        
                % all step errors
                errs_allsteps_gauss =  [errs_allsteps_gauss, errs_gauss]; % gauss
                % all frame errors
                errs_allframes_gauss = [errs_allframes_gauss,errs_gauss(end)]; % gauss
            %
            % plots    
                % gausssian
                f_gaussian = figure(22);  set(f_gaussian,'visible','off')
                subplot(121)
                imagesc(gaussian_dI); axis image, colormap(gray)
                title(sprintf('gaussian - frame %d/%d',i,length(multiFrame)))
                drawnow          
                subplot(122)
                plot(1:length(errs_allframes_gauss),errs_allframes_gauss,'Color',blu);
                xlabel('#frames','FontSize', 15), ylabel('residual |Fx - y| / pixel','FontSize', 15)
                title(sprintf('frame %d/%d',i,length(multiFrame)))
                drawnow    
                % save figure
                figure(f_gaussian); set(f_gaussian,'visible','off')
                filename = ['mfd_deconv_gaussian_with_curve'];
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)
                %----- gaussian deconved image -----
                f_gaussian = figure; set(f_gaussian,'visible','off');
                imagesc(clip(gaussian_dI,1,0)); axis image, colormap(gray)
                title('gaussian')
                filename = sprintf('deconv_gaussian_%d_frame',i);
                filename = fullfile(figPath,filename);
                print(gcf, '-depsc2', filename)
            end
            I = gaussian_dI;
        otherwise
            print '[mfd.m]: option.method can be one of those: 'pncg', 'cg', 'rl', 'gaussian''
    end
end
%##### figure #####
% for debug
figure(fclk); set(fclk,'visible','on'),
subplot(121)
title('each frame');
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex','FontSize', 15)
xlabel('$time/sec$','Interpreter','Latex','FontSize', 15)
subplot(122)
title('each frame');
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex','FontSize', 15); 
xlabel('$time/sec$','Interpreter','Latex','FontSize', 15)
figure(fstp); set(fstp,'visible','on'),
subplot(121)
title('each frame');
% legend('my pncg','my cg','my lucy','my gaussian');legend boxoff
ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex','FontSize', 15)
xlabel('$\#steps$','Interpreter','Latex','FontSize', 15)
subplot(122)
title('each frame');
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex','FontSize', 15);
xlabel('$\#steps$','Interpreter','Latex','FontSize', 15)
% for latex
saveResultFigure;
%%
% figure ||F*x - y||

fmfd_allStp = figure; clf,set(fmfd_allStp,'visible','off'),
loglog(errs_allsteps_pncg,'Color',dre),hold on
loglog(errs_allsteps_cg,'Color',mpg),hold on
loglog(errs_allsteps_lucy,'Color',ora),hold on
loglog(errs_allsteps_gauss,'Color',blu),hold on
xlabel('$\#steps$','Interpreter','Latex','FontSize', 15),ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex','FontSize', 15)
% legend('pncg','cg','lucy','gauss')
title('multi-frame')
filename = ['compAllStp_mfd'];
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

end