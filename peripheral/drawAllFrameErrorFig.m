function drawAllFrameErrorFig(residualErrAllFrame, relativeErrAllFrame, numFrame, counter, optimizerName, color, figPath, mode, fhandel)
%
% -------- ground truth frame error figure --------
switch mode
    case 'latex'
        % for latex
        % residual error
        figure;  set(gcf,'visible','off'),
        hData   = plot(0:length(residualErrAllFrame)-1,residualErrAllFrame,'Color', color);          
        hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
        hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        hTitle  = title(sprintf('%d frames', numFrame));
        set(gca,'XTick',0:length(residualErrAllFrame)-1)
        thisFigure;   
        drawnow
        filename = sprintf(strcat('mbd_residualErrorAllFrame_',optimizerName));
        filename = fullfile(figPath,filename);
        print(gcf, '-depsc2', filename)
        close gcf;
        % relative error
        figure;  set(gcf,'visible','off'),
        hData = plot(0:length(relativeErrAllFrame)-1,relativeErrAllFrame,'Color', color);          
        hYLabel = ylabel('$relative\ error$', 'Interpreter','Latex');
        hXLabel = xlabel('$\#frames$', 'Interpreter','Latex');
        hTitle = title(sprintf('%d frames', numFrame));
        set(gca,'XTick',0:length(relativeErrAllFrame)-1)
        thisFigure;   
        drawnow
        filename = sprintf(strcat('mbd_relativeErrorAllFrame_',optimizerName));
        filename = fullfile(figPath,filename);
        print(gcf, '-depsc2', filename)
        close gcf;
    case 'debug'
        % for debug
        figure(fhandel); subplot(121)
        hData = plot(0:length(residualErrAllFrame)-1,residualErrAllFrame,'Color', color);          
        hYLabel = ylabel('$\|Fx - y\| / pixel$', 'Interpreter','Latex');
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