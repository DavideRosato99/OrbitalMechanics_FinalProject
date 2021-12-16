function PorkchopPlot1( tPlot, costTot, costMin, tolCostMaxPlot, type )
% PorkchopPlot1 function that represent the porkchopplot
%
% PROTOTYPE:
%   PorkchopPlot1( tPlot, costTot, costMin, tolCostMaxPlot, type )
%
% INPUT:
%   tPlot[struct]       struct of dates
%   tPlot.Dep[1xN]      vector of departure dates [mjd2000]
%   tPlot.Arr[1xN]      vector of arrival dates [mjd2000]
%   costTot[NxN]        matrix of DV [km/s]
%   costMin[1]          minimum DV [km/s]
%   tolCostMaxPlot[1]   range of DV in the plot [km/s]
%   type[1]             representation option [ 'Dates' || 'ToF' ] [string]
%
% CONTRIBUTORS:
%   Fabio Spada
%
% VERSIONS:
%   2021-02-11
%

F = figure;
% F.Colormap = 'jet';

tDepPlot = tPlot.Dep;
tArrPlot = tPlot.Arr;
tArrMinCostPlot = tPlot.tArrMinCost;
tDepMinCostPlot = tPlot.tDepMinCost;
for ii = 1 : size(tDepPlot,2)
    tDepPlotn(ii) = datenum(mjd20002date(tDepPlot(ii)));
end
for ii = 1 : size(tArrPlot,2)
    tArrPlotn(ii) = datenum(mjd20002date(tArrPlot(ii)));
end

costTot(costTot > floor(costMin)+tolCostMaxPlot) = NaN;

if strcmp(type, 'Dates')
    
    [C, h] = contourf(tDepPlot, tArrPlot, costTot, floor(costMin) + linspace(0, tolCostMaxPlot-1+eps, tolCostMaxPlot) , 'showtext', 'off' );
    caxis(floor(min(costMin)) + [0 tolCostMaxPlot]);
    xtickangle(45);
    ytickangle(45);
    datetick('x', 'yyyy mmm dd', 'keeplimits');
    datetick('y', 'yyyy mmm dd', 'keeplimits');
    hold on
    set(gca, 'TickLabelInterpreter', 'latex');
    
    hcb = colorbar;
    set(hcb, 'TickLabelInterpreter', 'latex');
    
    gca;
    [C2, h2] = contour(tDepPlot, tArrPlot, tArrPlot' - tDepPlot, [0 : max(max(tArrPlot'-tDepPlot))/5 : max(max(tArrPlot'-tDepPlot))], 'k'); % remember to use same inputs as contour plotting
    clabel(C2, h2);
    
    if ~isempty(tArrMinCostPlot)
%         plot(tDepMinCostPlot, tArrMinCostPlot, 'o', 'MarkerFaceColor', [255, 0, 0]/256);
    end
    
elseif strcmp(type, 'ToF')
    
    [tDepPlotn, tArrPlotn] = meshgrid(tDepPlotn, tArrPlotn);
%     [C, h] = contourf(tDepPlot, tArrPlot-tDepPlot, costTot, floor(costMin) + linspace(0, tolCostMaxPlot-1+eps, tolCostMaxPlot) , 'showtext', 'off' );
    sc = surfc(tDepPlotn, tArrPlotn-tDepPlotn, costTot,'EdgeColor','none');
%     S.EdgeAlpha = 0; 
%     S.FaceAlpha = 0.5;
    caxis(floor(min(costMin)) + [0 tolCostMaxPlot]);
    xtickangle(-45);
    ytickangle(45);
    datetick('x', 'yyyy mmm dd', 'keeplimits');
    ylim([0 max(max(tArrPlotn-tDepPlotn))])
    hold on
    Ax = gca;
    Ax.TickLabelInterpreter = 'latex';
    Ax.Colormap = jet(256);
    Ax.Color = 'none';
%     set(gca, 'TickLabelInterpreter', 'latex', 'Colormap','jet', 'Color', 'None');
    
    hcb = colorbar;
    set(hcb, 'TickLabelInterpreter', 'latex');
    
    if ~isempty(tArrMinCostPlot)
%         plot(tDepMinCostPlot, tArrMinCostPlot, 'o', 'MarkerFaceColor', [255, 0, 0]/256);
    end
    
end
end