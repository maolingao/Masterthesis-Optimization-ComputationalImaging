%%
% should be defined separately for each figure
%
%{
% Draw figures
hData1 = plot(x,y);
hData2 = plot(x,y);
hData3 = plot(x,y);
hData4 = plot(x,y);
color = dre;
% Add labels
hTitle = title('My Publication-Quality Graphics');
hXLabel = xlabel('Length (m)');
hYLabel = ylabel('Mass (kg)');
% Add legend
hLegend = legend([hData1, hData2, hData3, hData4, hData5], ...
    'Data ({\it\mu} \pm {\it\sigma})', 'Fit (C{\itx}^3)', ...
    'Validation Data', 'Model (C{\itx}^3)', '95% CI', );
% Add text
hText = text(10, 800, ...
    sprintf('{\\itC = %0.1g \\pm %0.1g (CI)}', c, cint(2)-c));
%}
%% basic setup;
% mpg = [0,0.4717,0.4604]; % color [0,125,122]
% dre = [0.4906,0,0]; % color [130,0,0]
% ora = [255,153,51] ./ 255;
% blu = [0,0,0.509];
% gra = 0.5 * ones(3,1);

%%
if exist('hData','var')
% if ~exist('color','var') && exist('hData','var')
    hDataInfo = get(hData);
    color = hDataInfo.Color;
    linestyle = hDataInfo.LineStyle;
end

if ~exist('location','var') && exist('hLegend','var')
    hLegendInfo = get(hLegend);
    location = hLegendInfo.Location;
end
%
% set(gcf,'visible','on');
if exist('hData', 'var')    
    % Adjust line properties (functional & aesthetics)
    set(hData, 'LineWidth', 1, 'LineStyle',linestyle, 'Color', color)
    % set(hData1, 'LineStyle', '-.', 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
    % set(hData2, 'LineStyle', ':' , 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
    % set(hData3, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
    % set(hData4, 'LineStyle', '-' , 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
end
if exist('hTitle', 'var')    
    % Adjust title properties (aesthetics)
    set(hTitle, 'FontName', 'AvantGarde', 'FontSize', 17 )    
end
if exist('hXLabel', 'var')   
    % Adjust xlabel properties (aesthetics) 
    set(hXLabel, 'FontName', 'AvantGarde', 'Interpreter','Latex', 'FontSize', 15 )    
end
if exist('hYLabel', 'var')  
    % Adjust ylabel properties (aesthetics)  
    set(hYLabel, 'FontName', 'AvantGarde', 'Interpreter','Latex', 'FontSize', 15 )    
end
if exist('hLegend', 'var')  
    % Adjust legend properties (aesthetics)  
    set(hLegend, 'FontName', 'Helvetica', 'color','none', 'Location', location, 'FontSize', 13 )    
end

% Adjust axes properties
set(gca, 'FontName', 'Helvetica','FontSize', 8 )
axis tight;
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
%%
clear hData hYLabel hXlabel hTitle hLegend
%% additional options
% Adjust line properties (aesthetics)
%{
set(hFit, 'LineWidth', 2)
set(hE, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
    'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7])
set(hData, 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
set(hModel, 'LineWidth', 1.5)
set(hCI(1), 'LineWidth', 1.5)
set(hCI(2), 'LineWidth', 1.5)
%}
