%%
% should be defined separately for each figure
%
%%
if exist('hData','var') && ~isempty(hData)
% if ~exist('color','var') && exist('hData','var')
hDataInfo = get(hData);
% color = hDataInfo.Color;
linestyle = hDataInfo.LineStyle;
end
if ~exist('location','var') && exist('hLegend','var')
hLegendInfo = get(hLegend);
location = hLegendInfo.Location;
end
%
% set(gcf,'visible','on');
if exist('hData', 'var') && ~isempty(hData)
% Adjust line properties (functional & aesthetics)
set(hData, 'LineWidth', 2, 'LineStyle',linestyle);%, 'Color', color)
% set(hData1, 'LineStyle', '-.', 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
% set(hData2, 'LineStyle', ':' , 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
% set(hData3, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
% set(hData4, 'LineStyle', '-' , 'Marker', '.', 'LineWidth', 2, 'LineStyle',linestyle, 'Color', color)
end
if exist('hTitle', 'var')
% Adjust title properties (aesthetics)
% set(hTitle, 'FontName', 'AvantGarde', 'FontSize', 17 )
end
if exist('hXLabel', 'var')
% Adjust xlabel properties (aesthetics)
% set(hXLabel, 'FontName', 'AvantGarde', 'Interpreter','Latex', 'FontSize', 15 )
set(hXLabel, 'Interpreter','Latex' )
end
if exist('hYLabel', 'var')
% Adjust ylabel properties (aesthetics)
% set(hYLabel, 'FontName', 'AvantGarde', 'Interpreter','Latex', 'FontSize', 15 )
set(hYLabel, 'Interpreter','Latex')
end
if exist('hLegend', 'var')
% Adjust legend properties (aesthetics)
% set(hLegend, 'FontName', 'Helvetica', 'color','none', 'Location', location, 'FontSize', 13 )
set(hLegend, 'color','none', 'Location', location)
end
% Adjust axes properties
set(gca, 'FontName', 'Helvetica','FontSize', 8 )
set(gca, 'Box', 'off', 'TickDir', 'in', 'TickLength', [.02 .02], ...
'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
'LineWidth', 1.8)
axis tight; box on;
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
