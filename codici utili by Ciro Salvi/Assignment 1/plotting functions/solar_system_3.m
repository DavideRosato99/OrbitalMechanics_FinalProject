function [s_sys_3] = solar_system_3( p1, p2, p3, R_1, R_2, R_3, T_t_1, T_t_2, T_t_3, T_t_4, r_t_1, r_t_2, r_t_3, r_t_4, v_p1, v_p21, v_p22, v_p3, Sun_fct, p1_fct, p2_fct, p3_fct, SOI_fct, tr_dotted ) 
%solar_system_3 function that represent three planet and the transfer orbit for the
% interplanetary mission
%
% PROTOTYPE:
%   [s_sys_3] = solar_system_3( p1, p2, p3, R_1, R_2, R_3, T_t_1, T_t_2, T_t_3, T_t_4, r_t_1, r_t_2, r_t_3, r_t_4, v_p1, v_p21, v_p22, v_p3, Sun_fct, p1_fct, p2_fct, p3_fct, SOI_fct, tr_dotted )
%
% INPUT:
%   p1 , p2, p3[1]  indices of the planets to be plotted
%                        1:   Mercury
%                        2:   Venus
%                        3:   Earth
%                        4:   Mars
%                        5:   Jupiter
%                        6:   Saturn
%                        7:   Uranus
%                        8:   Neptune
%                        9:   Pluto
%   R_1[1]            Mean distance from the Sun to the planet_1 [AU]
%   R_2[1]            Mean distance from the Sun to the planet_2 [AU]
%   R_3[1]            Mean distance from the Sun to the planet_3 [AU]
%   T_t_1[1]          Date of departure 1 [Julian date]
%   T_t_2[1]          Date of arrival 1 [Julian date] 
%   T_t_3[1]          Date of departure 2 [Julian date]
%   T_t_4[1]          Date of arrival 2 [Julian date] 
%   r_t_1[3x1]        starting position of the transfer trajectory 1 [km]
%   r_t_2[3x1]        ending position of the transfer trajectory 1 [km]
%   r_t_3[3x1]        starting position of the transfer trajectory 2 [km]
%   r_t_4[3x1]        ending position of the transfer trajectory 2 [km]
%   v_p1[3x1]         velocity of planet_1 at starting position of the transfer trajectory 1 [km/s] 
%   v_p21[3x1]        velocity of planet_2 at ending position of the transfer trajectory 1 [km/s] 
%   v_p22[3x1]        velocity of planet_2 at starting position of the transfer trajectory 2 [km/s]
%   v_p3[3x1]         velocity of planet_3 at ending position of thetransfer trajectory 2 [km/s]
%   Sun_fct[1]        size factor for Sun
%   p1_fct[1]         size factor for planet_1
%   p2_fct[1]         size factor for planet_2
%   p3_fct[1]         size factor for planet_3
%   SOI_fct[1]        size factor of SOI
%   tr_dotted[1]      [ 0 1 ] switch for a complete dotted transfer orbit
%
% OUTPUT:
%   r_SOI_1[1]        radius of SOI of the planet_1 [km]
%   r_SOI_2[1]        radius of SOI of the planet_2 [km] 
%   r_SOI_3[1]        radius of SOI of the planet_3 [km] 
%   e_t1[1]           eccentricity of the first transfer orbit [-] 
%   e_t2[1]           eccentricity of the second transfer orbit [-] 
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

% Physical parameter 
AU = 1.495978707e8;                % [km] Astronomical Unit
% Sun
mu_S = astroConstants(4);          % [km^3/s^2] Sun's planetary constants
mean_radius_S = astroConstants(3); % [km] mean radius Sun
% Planet_1
mu_p1 = astroConstants(p1+10);     % [km^3/s^2] Planet_1's planetary constants
mr_p1 = astroConstants(p1+20);     % [km] mean radius Planet_1
[kep_1_arr(:),~] = uplanet( T_t_1, p1 );
a_p1 = kep_1_arr(1);               % [km] semimajor axis Planet_1
% Planet_2
mu_p2 = astroConstants(p2+10);     % [km^3/s^2] Planet_2's planetary constants
mr_p2 = astroConstants(p2+20);     % [km] mean radius Planet_2
[kep_2_arr(:),~] = uplanet( T_t_1, p2 );
a_p2 = kep_2_arr(1);               % [km] semimajor axis Planet_2
% Planet_3
mu_p3 = astroConstants(p3+10);     % [km^3/s^2] Planet_3's planetary constants
mr_p3 = astroConstants(p3+20);     % [km] mean radius Planet_3
[kep_3_arr(:),~] = uplanet( T_t_1, p3 );
a_p3 = kep_3_arr(1);               % [km] semimajor axis Planet_3

Ratio_s1 = mean_radius_S / mr_p1;  % ratio between Planet_1's radius and Planet_2's radius
Ratio_12 = mr_p1 / mr_p2;          % ratio between Planet_1's radius and Planet_2's radius
Ratio_13 = mr_p1 / mr_p3;          % ratio between Planet_1's radius and Planet_3's radius

r_SOI_1 = norm(R_1*AU)*(mu_p1/mu_S)^(2/5); % [km] radius of the Sphere of Influence
r_SOI_2 = norm(R_2*AU)*(mu_p2/mu_S)^(2/5); % [km] radius of the Sphere of Influence
r_SOI_3 = norm(R_3*AU)*(mu_p3/mu_S)^(2/5); % [km] radius of the Sphere of Influence

%% Orbit 
% Set options
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

if( nargin < 24 )
   tr_dotted = 0;end

% Set transfer trajectory 1 data
Dt1 = T_t_2 - T_t_1;                  % [days] duration of the transfer orbit 1
tspan_t1 = [ 0 Dt1*3600*24 ];         % arc of ellipse of transfer orbit 1
[a_t1,~,e_t1,~,v_t_1,v_t_2,~,~] = lambertMR( r_t_1, r_t_2, Dt1*24*3600 , mu_S, 0, 0, 0, 1);
if( e_t1<1 && e_t1>0 && tr_dotted==1 )
    T_t1 = (2*pi*sqrt(a_t1^3/mu_S));      % [sec] period of the transfer orbit 1
    tspan_t_dot1 = [ 0 T_t1-tspan_t1(2)]; % dotted arc of ellipse 1
end
y0_t1     = [ r_t_1; v_t_1' ];        % Initial conditions of the transfer orbit 1
y0_t_dot1 = [ r_t_2; v_t_2' ];        % Final conditions of the transfer orbit 1

% Set transfer trajectory 2 data
Dt2 = T_t_4 - T_t_3;                  % [days] duration of the transfer orbit 2
tspan_t2 = [ 0 Dt2*3600*24 ];         % arc of ellipse of transfer orbit 2
[a_t2,~,e_t2,~,v_t_3,v_t_4,~,~] = lambertMR( r_t_3, r_t_4, Dt2*24*3600 , mu_S, 0, 0, 0, 1);
if( e_t2<1 && e_t2>0 && tr_dotted==1 )
    T_t2 = (2*pi*sqrt(a_t2^3/mu_S));      % [sec] period of the transfer orbit 2
    tspan_t_dot2 = [ 0 T_t2-tspan_t2(2)]; % dotted arc of ellipse 2
end
y0_t2     = [ r_t_3; v_t_3' ];        % Initial conditions of the transfer orbit 2
y0_t_dot2 = [ r_t_4; v_t_4' ];        % Final conditions of the transfer orbit 2

% Set planets data
T_p1 = (2*pi*sqrt(a_p1^3/mu_S));  % [sec] period of Planet_1
T_p2 = (2*pi*sqrt(a_p2^3/mu_S));  % [sec] period of Planet_2
T_p3 = (2*pi*sqrt(a_p3^3/mu_S));  % [sec] period of Planet_2
tspan_p1 = [ 0 T_p1 ];
tspan_p2 = [ 0 T_p2 ];
tspan_p3 = [ 0 T_p3 ];
y0_p1 = [ r_t_1; v_p1 ];          % Initial conditions of Planet_1
y0_p2 = [ r_t_2; v_p22 ];         % Initial conditions of Planet_2
y0_p3 = [ r_t_4; v_p3 ];          % Initial conditions of Planet_3

% Perturbation of the orbits
[~, Y_1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p1, y0_p1, options);
[~, Y_2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p2, y0_p2, options);
[~, Y_3] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p3, y0_p3, options);
[~, Y_t1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t1, y0_t1, options);
if( e_t1<1 && e_t1>0 && tr_dotted==1 )
[~, Y_t_dot1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t_dot1, y0_t_dot1, options);end
[~, Y_t2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t2, y0_t2, options);
if( e_t2<1 && e_t2>0 && tr_dotted==1 )
[~, Y_t_dot2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t_dot2, y0_t_dot2, options);end

%evaluate where is the Planet 1 at arrival time 1
[kep_1_arr(:),~] = uplanet( T_t_2, p1 );
    [ r_1_arr1(:),~ ] = kep2car( kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6), mu_S );
%evaluate where is the Planet 2 at departure time 1
[kep_2_dep(:),~] = uplanet( T_t_1, p2 );
    [ r_2_dep1(:),~] = kep2car( kep_2_dep(1), kep_2_dep(2), kep_2_dep(3), kep_2_dep(4), kep_2_dep(5), kep_2_dep(6), mu_S );
%evaluate where is the Planet 3 at departure time 1
[kep_3_dep(:),~] = uplanet( T_t_1, p3 );
    [ r_3_dep1(:),~] = kep2car( kep_3_dep(1), kep_3_dep(2), kep_3_dep(3), kep_3_dep(4), kep_3_dep(5), kep_3_dep(6), mu_S );
%evaluate where is the Planet 3 at arrival time 1
[kep_3_arr(:),~] = uplanet( T_t_2, p3 );
    [ r_3_arr1(:),~] = kep2car( kep_3_arr(1), kep_3_arr(2), kep_3_arr(3), kep_3_arr(4), kep_3_arr(5), kep_3_arr(6), mu_S );

%evaluate where is the Planet 1 at departure time 2
[kep_1_dep(:),~] = uplanet( T_t_3, p1 );
    [ r_1_dep2(:),~] = kep2car( kep_1_dep(1), kep_1_dep(2), kep_1_dep(3), kep_1_dep(4), kep_1_dep(5), kep_1_dep(6), mu_S );
%evaluate where is the Planet 1 at arrival time 2
[kep_1_arr(:),~] = uplanet( T_t_4, p1 );
    [ r_1_arr2(:),~] = kep2car( kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6), mu_S );
%evaluate where is the Planet 2 at arrival time 2
[kep_2_arr(:),~] = uplanet( T_t_4, p2 );
    [ r_2_arr2(:),~] = kep2car( kep_2_arr(1), kep_2_arr(2), kep_2_arr(3), kep_2_arr(4), kep_2_arr(5), kep_2_arr(6), mu_S );
%evaluate where is the Planet 3 at departure time 2
[kep_3_dep(:),~] = uplanet( T_t_3, p3 );
    [ r_3_dep2(:),~] = kep2car( kep_3_dep(1), kep_3_dep(2), kep_3_dep(3), kep_3_dep(4), kep_3_dep(5), kep_3_dep(6), mu_S );
%% Plot Orbit
if( nargin < 19 ) % with no size factors the plot will be correctly scaled
    Sun_fct = 1;
    p1_fct  = 1;
    p2_fct  = 1;
    p3_fct  = 1;
end
% figure;
% figure('WindowState','maximized','Color','black'); %'WindowState','maximized',
% [0.07 0.62 1] [1.00 0.41 0.16] [0.39 0.83 0.07] [1.00,0.82,0.41]
hold on; grid on
plot3 ( Y_1(:,1),Y_1(:,2),Y_1(:,3), 'Color',[0.07 0.62 1],'LineWidth',1);    % Planet_1 orbit '#0072BD',
plot3 ( Y_2(:,1),Y_2(:,2),Y_2(:,3), 'Color',[1.00 0.41 0.16],'LineWidth',1);    % Planet_2 orbit '#D95319',
plot3 ( Y_3(:,1),Y_3(:,2),Y_3(:,3), 'Color',[0.39 0.83 0.07],'LineWidth',1);    % Planet_3 orbit '#4DBEEE'
if( e_t1<1 && e_t1>0 && tr_dotted==1 )
plot3 ( Y_t_dot1(:,1), Y_t_dot1(:,2), Y_t_dot1(:,3), 'Color','#77AC30','LineStyle','--','LineWidth',1);end % Transfer orbit 1 dotted
plot3 ( Y_t1(:,1),Y_t1(:,2),Y_t1(:,3), 'Color','#4DBEEE','LineWidth',1); % Transfer orbit 1
if( e_t2<1 && e_t2>0 && tr_dotted==1 )
plot3 ( Y_t_dot2(:,1), Y_t_dot2(:,2), Y_t_dot2(:,3), 'Color','#EDB120','LineStyle','--','LineWidth',1);end % Transfer orbit 2 dotted
plot3 ( Y_t2(:,1),Y_t2(:,2),Y_t2(:,3), 'Color',[1.00,0.82,0.41],'LineWidth',1); % Transfer orbit 2 '#EDB120'

factor_p1  = p1_fct * Sun_fct/Ratio_s1;                    % Size factor for planet_1
Planet_plot ( 10,[ 0 0 0 ], Sun_fct );                     % Sun plot

Planet_plot (  p1, r_t_1,    factor_p1 );                  % Planet_1 plot departure 1
Planet_plot (  p1, r_1_arr1, factor_p1 );                  % Planet_1 plot arrival 1
Planet_plot (  p1, r_1_dep2, factor_p1 );                  % Planet_1 plot departure 2
Planet_plot (  p1, r_1_arr2, factor_p1 );                  % Planet_1 plot arrival 2

Planet_plot (  p2, r_2_dep1, p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot departure 1
Planet_plot (  p2, r_t_2,    p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot arrival 1
Planet_plot (  p2, r_t_3,    p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot departure 2
Planet_plot (  p2, r_2_arr2, p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot arrival 2

Planet_plot (  p3, r_3_dep1, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot departure 1
Planet_plot (  p3, r_3_arr1, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot arrival 1
Planet_plot (  p3, r_3_dep2, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot departure 2
Planet_plot (  p3, r_t_4,    p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot arrival 2
set(gca,'color','none','TickLabelInterpreter','latex') % no background for the graph

view(35,20);
axis equal
if( nargin >= 23 ) % SOI plot
SOI_plot ( r_SOI_2*SOI_fct*Sun_fct, 0.5, r_t_2);end


s_sys_3.r_SOI_1 = r_SOI_1;
s_sys_3.r_SOI_2 = r_SOI_2;
s_sys_3.r_SOI_3 = r_SOI_3;
s_sys_3.e_t1 = e_t1;
s_sys_3.e_t2 = e_t2;
end