% Interplanetary_mission starting from the overall window this script
% represent the porkchop plot for each leg

close all; clear all; clc

% Physical parameter 
mu_S = astroConstants(4);   % [km^3/s^2] Sun's planetary constants
mu_M = astroConstants(14);  % [km^3/s^2] Marte's planetary constants
mu_V = astroConstants(13);  % [km^3/s^2] Earth's planetary constants
mu_J = astroConstants(12);  % [km^3/s^2] Venus's planetary constants
AU = astroConstants(2);         % [km] Astronomical Unit
% r_M = 0.4;                  % [AU] Distance from the Sun
% r_V = 0.7;                  % [AU] Distance from the Sun
% r_J = 5.2;                  % [AU] Distance from the Sun

% Choice of the planet
% 4 Marte - 3 Earth - 2 venus
% p1 = 1;                     % Mercury
p1 = 4;                     % Venus
p2 = 3;                    % Jupiter

max_dv = 50;

mean_radius_S = astroConstants(3);
mean_radius_M = astroConstants(21); r_M = mean_radius_M;
mean_radius_E = astroConstants(22);
mean_radius_V = astroConstants(25);
Ratio_rs = mean_radius_S / mean_radius_M; % ratio between Sun's radius and Planet_1's radius
Ratio_r = mean_radius_M / mean_radius_E; % ratio between Planet_1's radius and Planet_2's radius
% Ratio_r = mean_radius_V / mean_radius_J; % ratio between Planet_1's radius and Planet_2's radius

t_start_1 = [ 2027, 12, 1, 0, 0, 0];
t_end_1   = [ 2042, 12, 1, 0, 0, 0];
t_start_2 = [ 2027, 12, 1, 0, 0, 0];
t_end_2   = [ 2042, 12, 1, 0, 0, 0];

% Transform the date in Modified Julian day 2000 number from Gregorian date
t_start_j_1 = date2mjd2000(t_start_1);
t_end_j_1   = date2mjd2000(t_end_1);
t_start_j_2 = date2mjd2000(t_start_2);
t_end_j_2   = date2mjd2000(t_end_2);
% reverse fun: mjd20002date(a)

% T_1 = linspace( t_start_j_1, t_end_j_1, (t_end_j_1 - t_start_j_1 + 1)); % Date range on the departure orbit
% T_2 = linspace( t_start_j_2, t_end_j_2, (t_end_j_2 - t_start_j_2 + 1)); % Date range on the arrival orbit
T_1 = linspace( t_start_j_1, t_end_j_1, 100); % Date range on the departure orbit
T_2 = linspace( t_start_j_2, t_end_j_2, 100); % Date range on the arrival orbit

%% Evaluate orbital parameters
for ii = 1 : length(T_1)
    [kep_1(ii,:)] = uplanet( T_1(ii), p1 );
    orb_1 = [kep_1(ii,1), kep_1(ii,2), kep_1(ii,3), kep_1(ii,4), kep_1(ii,5), kep_1(ii,6)];
    [ r_1(ii,:),v_1(ii,:) ] = par2car(orb_1 , mu_S );
end
for ii = 1 : length(T_2)
    [kep_2(ii,:)] = uplanet( T_2(ii), p2 );
    orb_2 = [kep_2(ii,1), kep_2(ii,2), kep_2(ii,3), kep_2(ii,4), kep_2(ii,5), kep_2(ii,6)];
    [ r_2(ii,:),v_2(ii,:) ] = par2car( orb_2, mu_S );
end
DV1 = zeros( length(T_1), length(T_2) );
DV2 = zeros( length(T_1), length(T_2) );

z = 0;
for ii = 1 : length(T_1)
    for jj = 1 : length(T_2)
        z = z +1 ;
        if( T_2(jj) > T_1(ii) )            % check if arrival date is more than departure date
            Dt(ii,jj) = T_2(jj) - T_1(ii); % [days]
            [~,~,~,~,v_1_l,v_2_l] = lambertMR( r_1(ii,:), r_2(jj,:), Dt(ii,jj)*24*3600 , mu_S, 0, 0, 0, 1);
            DV1(ii,jj) = norm( v_1_l - v_1(ii,:) );
            DV2(ii,jj) = norm( v_2(jj,:) - v_2_l );
			V_1_L(ii,jj,:) = v_1_l;   
            V_2_L(ii,jj,:) = v_2_l;
            if ( DV1(ii,jj)>=max_dv )
                DV1(ii,jj) = 0;
                DV2(ii,jj) = 0;
            end
        end
    end
    fprintf('ii = %.2f%%\n', (ii*100/length(T_1)) ); % percentage of process    
end
DV = DV1 + DV2;
%% save data
out.DV_tot = DV;
out.DV_1   = DV1;
out.DV_2   = DV2;
out.T_1    = T_1;
out.T_2    = T_2;
out.v_p1   = v_1;
out.v_p2   = v_2;
out.V_L_1  = V_1_L;
out.V_L_2  = V_2_L;
beep
%% Transform the date in Gregorian date from Modified Julian day 2000 number 
for ii = 1 : length(T_1)
    T_1_g(ii) = datenum(mjd20002date(T_1(ii)));end
for ii = 1 : length(T_2)
    T_2_g(ii) = datenum(mjd20002date(T_2(ii)));end
%% Min value
[i,j] = find( DV1 == min(DV1(DV1>0)) );
% [i,j] = find( DV == min(DV(:)) );

t_departure_opt = datetime(mjd20002date(T_1(i)));
t_arrival_opt   = datetime(mjd20002date(T_2(j)));
fprintf('Minimum Delta_v = %.2f [km/s];\ndeparture date: %s\nArrival date: %s\n',DV(i,j),t_departure_opt,t_arrival_opt);
%% Plot
%% Porkchop
Porckchop_switch = 0;
if(Porckchop_switch)
figure; hold on; grid minor;
title('Porkchop');
contour( T_1_g , T_2_g , DV1', floor(min(min(DV1(DV1>0))))+(0:0.5:10), 'ShowText', 'off' );
caxis(floor(min(min(DV1(DV1>0))))+[0 15]);
caxis('manual');
xtickangle(45)
ytickangle(45)
datetick('x', 'yyyy mmm dd','keeplimits');
datetick('y', 'yyyy mmm dd','keeplimits');
ylabel('Arrival date');
xlabel('departure date');
hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
ha = gca;
ha.FontSize = 13;
hcb.Title.FontSize = 15;
% [C2,h2] = contour( T_1_g , T_2_g , T_2_g'-T_1_g, [60 120 180 240 300], 'k' ); %constant tof lines
% clabel(C2,h2);
end
%% Surf
surf_switch = 0;
if(surf_switch)
M = zeros(size(DV1,1), size(DV1,2));
for ii = 1 : size(DV1,1)
    for jj = 1 : size(DV1,2)
        if( DV1(ii,jj) < 15 )
            M(ii,jj) = DV1(ii,jj);
        end
    end
end
figure
mesh( T_1_g , T_2_g , M' );
datetick('x', 'yyyy mmm dd','keeplimits');
datetick('y', 'yyyy mmm dd','keeplimits');
ylabel('Arrival date');
xlabel('departure date');
colorbar
view(-15,25)
end
%% Orbit
% Set options
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
tspan_t = [ 0 Dt(i,j)*3600*24 ]; % arc of ellipse 
[A_t,~,~,~,v_1_l,v_2_l,~,~] = lambertMR( r_1(i,:), r_2(j,:), Dt(i,j)*24*3600 , mu_S, 0, 0, 0, 1);
tspan_t_dot = [ 0 ((2*pi*sqrt(abs(A_t)^3/mu_S))-(Dt(i,j)*3600*24))];
tspan_1 = [ 0 (2*pi*sqrt(kep_1(1,1)^3/mu_S))];
tspan_2 = [ 0 (2*pi*sqrt(kep_2(1,1)^3/mu_S))];
y0_t = [ r_1(i,:)'; v_1_l' ];
y0_t_dot = [ r_2(j,:)'; v_2_l' ];
y0_E = [ r_1(i,:)'; v_1(i,:)' ];
y0_M = [ r_2(j,:)'; v_2(j,:)' ];
[time_t, Y_t] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t, y0_t, options);
[time_t_dot, Y_t_dot] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t_dot, y0_t_dot, options);
[time_1, Y_1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_1, y0_E, options);
[time_2, Y_2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_2, y0_M, options);
%evaluate where is the Planet 1 at arrival time
[kep_1_arr(:),~] = uplanet( T_2(j), p1 );
orb_1_arr = [kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6)];
    [ r_1_arr(:),v_1_arr(:) ] = par2car( orb_1_arr, mu_S );
%evaluate where is the Planet 2 at departure time
[kep_2_arr(:),~] = uplanet( T_1(i), p2 );
orb_2_arr = [kep_2_arr(1), kep_2_arr(2), kep_2_arr(3), kep_2_arr(4), kep_2_arr(5), kep_2_arr(6)];
    [ r_2_dep(:),v_2_dep(:) ] = par2car( orb_2_arr, mu_S );
%% Plot Orbit Opt
Orbit_switch = 0;
if(Orbit_switch)
figure;hold on; grid on
% p_D = plot3 ( r_1(:,1), r_1(:,2), r_1(:,3), 'Color','#0072BD','LineWidth',15 ); %Departure window
% p_A = plot3 ( r_2(:,1), r_2(:,2), r_2(:,3), 'Color','#D95319','LineWidth',15 ); %Arrival window
% p_D.Color(4) = 0.3;
% p_A.Color(4) = 0.3;
plot3 ( Y_t_dot(:,1), Y_t_dot(:,2), Y_t_dot(:,3), 'Color','#77AC30','LineStyle','--','LineWidth',2);
plot3 ( Y_t(:,1),Y_t(:,2),Y_t(:,3), 'Color','#77AC30','LineWidth',2);
plot3 ( Y_1(:,1),Y_1(:,2),Y_1(:,3), 'Color','#0072BD','LineWidth',2);
plot3 ( Y_2(:,1),Y_2(:,2),Y_2(:,3), 'Color','#D95319','LineWidth',2);

factor_Sun = 40;                      % factor size for Sun
factor_p1  = 5e3*factor_Sun/Ratio_rs; % factor size for planet_1
Planet_plot ( 10,[ 0 0 0 ], factor_Sun );                                     % Sun plot
Planet_plot (  p2,[ r_2(j,1), r_2(j,2), r_2(j,3)], factor_p1/Ratio_r );       % Planet_2 plot arrival
Planet_plot (  p2,[ r_2_dep(1), r_2_dep(2), r_2_dep(3)], factor_p1/Ratio_r ); % Planet_2 plot departure
Planet_plot (  p1,[ r_1(i,1), r_1(i,2), r_1(i,3)], 3*factor_p1 );             % Planet_1 plot departure
Planet_plot (  p1,[ r_1_arr(1), r_1_arr(2), r_1_arr(3)], 3*factor_p1 );       % Planet_1 plot arrival
% legend({'Departure window','Arrival window'});
% view(2);
view(35,20);
axis equal
end
%% SOI
r_SOI = norm(r_M*AU)*(mu_V/mu_S)^(2/5); % [km] radius of the Sphere of Influence

check_SOI = 0;
if( check_SOI )
    figure; grid on; hold on;
    factor = 1;                        % factor on Earth-Sun distance
    Planet_plot ( 10,[ 0 0 0 ], 1 );  % Sun plot
    SOI_plot(r_SOI*10*factor, 0.5);    % SOI
    axis equal
    view(25,15)
end
p = profile ('info');