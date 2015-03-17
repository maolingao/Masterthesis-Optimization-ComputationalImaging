function drawAllFrameErrorFig(residualErrAllFrame, relativeErrAllFrame, numFrame, counter, optimizerName, color, figPath, mode, fhandel, kernelerror)
%
if ~exist('kernelerror','var')
    kernelerror = 0;
end
% -------- ground truth frame error figure --------
switch mode
    case 'latex'
        % for latex
        % residual error
        figure;  set(gcf,'visible','off'),
        hData   = plot(0:length(residualErrAllFrame)-1,residualErrAllFrame,'Color', color);          
        hYLabel = ylabel('$\|f * x - y\| / pixel$', 'Interpreter','Latex');
        hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        hTitle  = title(sprintf('%d frames', numFrame));
%         set(gca,'XTick',0:length(residualErrAllFrame)-1)
        if kernelerror == 1  % residual error of kernel of all frames
            filename = 'mbd_f_residualErrorAllFrame_';
            hXLabel = xlabel('$\#steps$', 'Interpreter','Latex');
        else
            filename = 'mbd_gt_residualErrorAllFrame_';
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        end
        thisFigure;  
%         if ~isempty(residualErrAllFrame)
%             ylim([0, max(residualErrAllFrame)]);
% %             xlim([0, 48]);
%         end
        drawnow
        filename = sprintf(strcat(filename,optimizerName));
        filename = fullfile(figPath,filename);
        print(gcf, '-depsc2', filename)
        close gcf;
        % relative error
        figure;  set(gcf,'visible','off'),
        hData = plot(0:length(relativeErrAllFrame)-1,relativeErrAllFrame,'Color', color);          
        hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
        hTitle = title(sprintf('%d frames', numFrame));
%         set(gca,'XTick',0:length(relativeErrAllFrame)-1)
        if kernelerror == 1  % relative error of kernel of all frames
            filename = 'mbd_f_relativeErrorAllFrame_';
            hXLabel = xlabel('$\#steps$', 'Interpreter','Latex');
        else
            filename = 'mbd_gt_relativeErrorAllFrame_';
            hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        end
        thisFigure;   
%         if ~isempty(residualErrAllFrame)
%             ylim([0, 1]);
% %             xlim([0, 48]);
%         end
        drawnow
        filename = sprintf(strcat(filename,optimizerName));
        filename = fullfile(figPath,filename);
        print(gcf, '-depsc2', filename)
        close gcf;
    case 'debug'
        % for debug
        figure(fhandel); subplot(121)
        hData = plot(0:length(residualErrAllFrame)-1,residualErrAllFrame,'Color', color);          
        hYLabel = ylabel('$\|f * x - y\| / pixel$', 'Interpreter','Latex');
        hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        hTitle = title(optimizerName);
        thisFigure;   
        drawnow
        subplot(122) 
        hData = plot(0:length(relativeErrAllFrame)-1,relativeErrAllFrame,'Color', color);       
        hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
        hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        hTitle = title(sprintf('frame %d/%d',counter,numFrame));
        thisFigure;   
        drawnow     
    otherwise
        disp('[drawAllFrameErrorFig.m] : mode can either be latex or debug.')
end
        
end