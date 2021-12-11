function [h1, h2] = porkChopPlot2DoF(dep_ID, arr_ID, dep_min, dep_max, arr_min, arr_max, dt_MJD, arcNum)
% PROTOTYPE:
% [h1, h2] = porkChopPlot2DoF(dep_ID, arr_ID, dep_min, dep_max, arr_min, arr_max, dt_MJD, arcNum)  
% 
% DESCRIPTION:
%   Returns one or two porkchop plots of one or two Lambert arcs
% INPUT:
%     dep_ID,arr_ID: Integer number identifying the celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%     dep_min, dep_max, arr_min, arr_max: dates vectors [year, month, day, 
%                       hours, minutes, seconds]
%     dt_MJD: time step-size
%     arcNum: number of Lambert arcs {1 or 2}
% 
% OUTPUT:
%     figure h1
%     figure h2

% CALLED FUNCTIONS:
%     astroConstants.m
%     date2mjd2000.m (time)
%     date2string.m
%     kep2car.m

mu_s = astroConstants(4); % [km^3/s^2]

sun_ID = 10;


launch_W = [dep_min:dt_MJD:dep_max];
arrival_W = [arr_min:dt_MJD:arr_max];

DeltaV_tot_mat = zeros(length(launch_W), length(arrival_W));
DeltaV_T1_mat = zeros(size(DeltaV_tot_mat));
DeltaV_T2_mat = zeros(size(DeltaV_tot_mat));


lastMsgNumel = 0; % Initialize counter of last message-to-delete elements (right before cycling)
for k = 1:length(launch_W)
    t1 = launch_W(k);
    
    kep_dep = uplanet(t1, dep_ID);
    [RR1, VV1] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_s);
    
    for j = 1:length(arrival_W)
        
        t2 = arrival_W(j);
        ToF_MJD2000 = t2-t1; % [Julian Days]
        
        if ToF_MJD2000 > 0
            ToF = ToF_MJD2000 * (24*3600); % [s]
            kep_arr = uplanet(t2, arr_ID);
            [RR2, VV2] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_s);

            % Lambert solver call
            [a,p,e,ERROR,VVT1,VVT2,TPAR,theta] = lambertMR( RR1, RR2 , ToF, mu_s, 0, 0, 2 );
            VVT1 = VVT1(:); VVT2 = VVT2(:);

            DeltaV_T1 = norm(VVT1-VV1);
            DeltaV_T2 = norm(VV2-VVT2);
            DeltaV_tot = DeltaV_T1 + DeltaV_T2;
            DeltaV_tot_mat(k,j) = DeltaV_tot;
            DeltaV_T1_mat(k,j) = DeltaV_T1;
            DeltaV_T2_mat(k,j) = DeltaV_T2;


            fprintf(repmat('\b',1,lastMsgNumel)); % delete latest progress-message
            progressString = sprintf('Progress: %.1f/100', floor(1000*((k-1)*length(arrival_W)+j)/(length(launch_W)*length(arrival_W)))/10);
            fprintf(progressString); % print new progress-message
            lastMsgNumel = numel(progressString); % update elements counter of last message-to-delete
        else
            DeltaV_T1_mat(k,j) = NaN;
            DeltaV_T2_mat(k,j) = NaN;
            DeltaV_tot_mat(k,j) = NaN;
        end
            
    end
end


    if arcNum == 1
        DeltaV_tot_mat = DeltaV_tot_mat - DeltaV_T2_mat;
    elseif arcNum == 2
        DeltaV_tot_mat = DeltaV_tot_mat - DeltaV_T1_mat;
    end

%% Pork-Chop Plot
h1 = figure;
noRepresThreshold = 100; % [km/s]
DeltaV_tot_mat(DeltaV_tot_mat > noRepresThreshold) = noRepresThreshold;

%subplot(1,2,1)
shift_MJD2000 = datenum(2000,01,01,12,00,00);
[AX1, AX2] = meshgrid(launch_W+shift_MJD2000, arrival_W+shift_MJD2000);
surf(AX1, AX2, DeltaV_tot_mat', 'LineStyle', 'None');
xlabel('Departure date')
ylabel('Arrival date')
datetick('x','yyyy mmm dd')
xlim(shift_MJD2000+[min(launch_W), max(launch_W)]);
datetick('y','yyyy mmm dd')
ylim(shift_MJD2000+[min(arrival_W), max(arrival_W)]);
zlabel('\DeltaV [km/s]')
hcb = colorbar;
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
shading interp
view(0,90)  % XY-plane view

%%
%subplot(1,2,2)
h2 = figure;
hold on
cLinesNum = 2000;
contour(AX1,AX2,DeltaV_tot_mat', cLinesNum, 'LineWidth', 1);
xlabel('Departure date')
ylabel('Arrival date')
datetick('x','yyyy mmm dd')
xlim(shift_MJD2000+[min(launch_W), max(launch_W)]);
datetick('y','yyyy mmm dd')
ylim(shift_MJD2000+[min(arrival_W), max(arrival_W)]);
axis equal
hcb = colorbar;
hcb.Limits = [min(min(DeltaV_tot_mat)), max(max(DeltaV_tot_mat))];
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
%myMap = colormap(parula);
%colormap(myMap);

% [contourMatrix,contourObject] = contour(AX1,AX2,AX2-AX1, 12, 'm');
% contourObject.LevelList=round(contourObject.LevelList, 0); % round level labels to 0 dec. digits
% clabel(contourMatrix,contourObject, 'Color','m', 'LabelSpacing',300, 'FontSize',11);


nLines = 12;
nLines = round(nLines/2)*2;
xAxis = AX1(1,:);
yAxis = AX2(:,1);
xRange = (xAxis(end)-xAxis(1));
yRange = (yAxis(end)-yAxis(1));
fSolveOptions = optimset('Display','off');
for k = -nLines/2:nLines/2
    y = @(x) yAxis(1) + 1*(x-xAxis(1)) + k*1.2*(yRange/xRange) * (yRange/nLines);
    y1 = @(x) y(x) - yAxis(1);
    yEnd = @(x) y(x) - yAxis(end);

    plot (xAxis, y(xAxis), 'm');
    x_at_y0 = fsolve(y,0,fSolveOptions);
    ToF = 0-x_at_y0;
    
    x_at_y1 = fsolve(y1,0,fSolveOptions);
    x_at_yEnd = fsolve(yEnd,0,fSolveOptions);

    % Determine ToF label position
    if x_at_y1 >= xAxis(1)
        x_textLab = (x_at_y1 + xAxis(end)) / 2;
    else
        x_textLab = (x_at_yEnd + xAxis(1)) / 2;
    end
    y_textLab = y(x_textLab);
    
    % Add ToF label only if it would be visible
    if x_textLab > xAxis(1) & x_textLab < xAxis(end) & y_textLab > yAxis(1) & y_textLab < yAxis(end);
        myText = text(x_textLab, y_textLab, num2str(round(ToF)));
        set(myText, 'Rotation',45, 'Color', 'm');
    end
    
end
 

 
 











end

