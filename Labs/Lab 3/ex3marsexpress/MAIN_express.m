
close all
clear all

addpath('time');    % Required for any call of ephNEO()
addpath('textures');


mu_sun = astroConstants(4); % [km^3/s^2]


dep_ID = 3; % Earth ID [-]
arr_ID = 4; % Mars ID [-]
sun_ID = 10;


% Launch Window
dep_min = [2003, 04, 01, 00, 00, 00];   % [year, month, day, hours, minutes, seconds]
dep_max = [2003, 08, 01, 23, 59, 59];

% Arrival Window
arr_min = [2003, 09, 01, 00, 00, 00];
arr_max = [2004, 03, 01, 23, 59, 59];


dep_min = date2mjd2000(dep_min);    dep_max = date2mjd2000(dep_max);
arr_min = date2mjd2000(arr_min);    arr_max = date2mjd2000(arr_max);

dt = date2mjd2000([2000, 01, 01, 20, 00, 00]); % 8 hours [MJD2000]


launch_W = [dep_min:dt:dep_max];
arrival_W = [arr_min:dt:arr_max];

DeltaV_tot_mat = zeros(length(launch_W), length(arrival_W));
DeltaV_T1_mat = zeros(size(DeltaV_tot_mat));
DeltaV_T2_mat = zeros(size(DeltaV_tot_mat));


v_inf_launcher = 2.969;
DeltaV_optConstr = inf;

lastMsgNumel = 0; % Initialize counter of last message-to-delete elements (right before cycling)
for k = 1:length(launch_W)
    t1 = launch_W(k);
    
    kep_dep = uplanet(t1, dep_ID);
    [RR1, VV1] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);
    
    for j = 1:length(arrival_W)
        
        t2 = arrival_W(j);
        ToF_MJD2000 = t2-t1; % [Julian Days]
        ToF = ToF_MJD2000 * (24*3600); % [s]
        kep_arr = uplanet(t2, arr_ID);
        [RR2, VV2] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);
        
        % Lambert solver call
        [a,p,e,ERROR,VVT1,VVT2,TPAR,theta] = lambertMR( RR1, RR2 , ToF, mu_sun, 0, 0, 2 );
        VVT1 = VVT1(:); VVT2 = VVT2(:);

        DeltaV_T1 = norm(VVT1-VV1);
        DeltaV_T2 = norm(VV2-VVT2);
        DeltaV_tot = DeltaV_T1 + DeltaV_T2;
        DeltaV_tot_mat(k,j) = DeltaV_tot;
        DeltaV_T1_mat(k,j) = DeltaV_T1;
        DeltaV_T2_mat(k,j) = DeltaV_T2;
        
        if k == 1 & j == 1
            % (i,j) = (1,1) -->   temporary assignment of optimum 
            %                     DeltaV value, for subsequent comparisons
            DeltaV_opt = DeltaV_tot;
            i_opt = k; j_opt = j; RR1_opt = RR1; RR2_opt = RR2;
        elseif DeltaV_tot < DeltaV_opt
                DeltaV_opt = DeltaV_tot;
                i_opt = k; j_opt = j; RR1_opt = RR1; RR2_opt = RR2;
        end
        
        
        if (DeltaV_T1 < v_inf_launcher) & (DeltaV_tot < DeltaV_optConstr)
            DeltaV_optConstr = DeltaV_tot;
            i_optConstr = k; j_optConstr = j; RR1_optConstr = RR1; RR2_optConstr = RR2;
        end
        
        fprintf(repmat('\b',1,lastMsgNumel)); % delete latest progress-message
        progressString = sprintf('Progress: %.1f/100', floor(1000*((k-1)*length(arrival_W)+j)/(length(launch_W)*length(arrival_W)))/10);
        fprintf(progressString); % print new progress-message
        lastMsgNumel = numel(progressString); % update elements counter of last message-to-delete
    end
end




%% Plot optimum Lambert's Arc



t1 = launch_W(i_opt);
t2 = arrival_W(j_opt);
ToF_MJD2000 = arrival_W(j_opt) - launch_W(i_opt); % [Julian Days]
ToF = ToF_MJD2000 * (24*3600); % [s]

kep_dep = uplanet(t1, dep_ID);
kep_arr = uplanet(t2, arr_ID);

[a,p,e,ERROR,VVT1,VVT2,TPAR,theta] = lambertMR( RR1_opt, RR2_opt , ToF, mu_sun, 0, 0, 2 );
RR1_opt = RR1_opt(:); RR2 = RR2(:); VVT1 = VVT1(:); VVT2 = VVT2(:);

y0 = [RR1_opt; VVT1];


fprintf('Minimum DeltaV arc w/ cost: %f km/s\n', DeltaV_tot_mat(i_opt, j_opt))
fprintf('Departure: %s\n', date2string(mjd20002date(t1)))
fprintf('Arrival: %s\n', date2string(mjd20002date(t2)))


% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

dt = 100; % Step size [s]

% Set time span
tspan = [0:dt:ToF];

% Perform the integration

% Lambert's arc
[time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X = StateMat_dep(:,1); Y = StateMat_dep(:,2); Z = StateMat_dep(:,3);
%VX = StateMat_dep(:,4); VY = StateMat_dep(:,5); VZ = StateMat_dep(:,6);

h = figure;
set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on
hold on
plot3(X,Y,Z, 'LineWidth', 3)



% Initial orbit
kep0 = uplanet(launch_W(1), dep_ID);
[rr0_dep, vv0_dep] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_sun);
y0 = [rr0_dep; vv0_dep];
T = 2*pi*sqrt(kep0(1)^3/mu_sun);
tspan = [0:dt:T];
[time_vec_dep, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X_dep = StateMat_dep(:,1); Y_dep = StateMat_dep(:,2); Z_dep = StateMat_dep(:,3);
%VX = StateMat_dep(:,4); VY = StateMat_dep(:,5); VZ = StateMat_dep(:,6);
plot3(X_dep,Y_dep,Z_dep, 'Color', [244, 66, 66]/255, 'LineStyle', '--')

% Final orbit
kep0 = uplanet(launch_W(1), arr_ID);
[rr0_arr, vv0_arr] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_sun);
y0 = [rr0_arr; vv0_arr];
T = 2*pi*sqrt(kep_arr(1)^3/mu_sun);
tspan = [0:dt:T];
[time_vec_arr, StateMat_arr] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X_arr = StateMat_arr(:,1); Y_arr = StateMat_arr(:,2); Z_arr = StateMat_arr(:,3);
%VX = StateMat_arr(:,4); VY = StateMat_arr(:,5); VZ = StateMat_arr(:,6);
plot3(X_arr,Y_arr,Z_arr, 'Color', [38, 180, 102]/255, 'LineStyle', '--')


view([198,27])
light('Position',[0 0 0],'Style','local');
sunGlobe = plotPlanet(sun_ID, [0,0,0], h, 30);
depGlobe = plotPlanet(dep_ID, RR1_opt, h, 20);
arrGlobe = plotPlanet(arr_ID, RR2_opt, h, 10);
axis equal

legend('Lambert''s Arc', 'Initial Orbit', 'Target Orbit')


%% Animated plot
h = figure;
hold on
grid on
plot3(X_dep,Y_dep,Z_dep, 'Color', [244, 66, 66]/255, 'LineStyle', '--')
plot3(X_arr,Y_arr,Z_arr, 'Color', [38, 180, 102]/255, 'LineStyle', '--')
view([198,27])
light('Position',[0 0 0],'Style','local');
sunGlobe = plotPlanet(sun_ID, [0,0,0], h, 30);
depDeltaTimeSeconds = (t1-launch_W(1))*24*3600;
depIndex = 0;
markerSize = 5;
j = 0; k = 0;
for k = 1:3000:length(time_vec_dep)
    
    if k > 1
       delete(depGlobe);    delete(arrGlobe);
    end
    
    depGlobe = plotPlanet(dep_ID, StateMat_dep(k,1:3), h, 20);
    arrGlobe = plotPlanet(arr_ID, StateMat_arr(k,1:3), h, 10);

    if time_vec_dep(k) >= depDeltaTimeSeconds
        if ~depIndex
            depIndex = k;
        end
        j = k-depIndex+1;
        if j > length(time_vec)
            break
        end
        
        if j > 1
            delete(SC);
            delete(traj_SC);
        end
        
        SC = plot3(X(j), Y(j), Z(j),'o',...    % Plot the current position of the S/C
                        'MarkerSize', markerSize,'MarkerEdgeColor','r','MarkerFaceColor',[0.8,0.2,0.2]);
        traj_SC = plot3(X(1:j), Y(1:j), Z(1:j), 'LineWidth', 1.5, 'Color', 'b');            
    end
    
    
    drawnow
    
        
    
end






%% Pork-Chop Plot
figure

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

%%
%subplot(1,2,2)
h = figure;
hold on
contour(AX1,AX2,DeltaV_tot_mat', 500, 'LineWidth', 1);
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
 

 
 








