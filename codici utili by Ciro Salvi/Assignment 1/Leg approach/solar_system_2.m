function [s_sys_2] = solar_system_2( p1, p2, R_1, R_2, T_t_1, T_t_2, r_t_1, r_t_2, v_p1, v_p2, Sun_fct, p1_fct, p2_fct, SOI_fct, tr_dotted ) 
%solar_system_2 function that represent two planet and the transfer orbit for the
% interplanetary mission
%
% PROTOTYPE:
%   [s_sys_2] = solar_system_2( p1, p2, R_1, R_2, T_t_1, T_t_2, r_t_1, r_t_2, v_p1, v_p2, Sun_fct, p1_fct, p2_fct, SOI_fct, tr_dotted )
%
% INPUT:
%   p1 , p2 [1]      indices of the planets to be plotted
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
%   T_t_1[1]          Date of departure [Julian date]
%   T_t_2[1]          Date of arrival [Julian date] 
%   r_t_1[3x1]        starting position of the transfer trajectory [km]
%   r_t_2[3x1]        ending position of the transfer trajectory [km]
%   v_p1[3x1]         velocity of planet_1 at starting position of the transfer trajectory [km/s] 
%   v_p2[3x1]         velocity of planet_2 at ending position of the transfer trajectory [km/s] 
%   Sun_fct[1]        size factor for Sun
%   p1_fct[1]         size factor for planet_1
%   p2_fct[1]         size factor for planet_2
%   SOI_fct[1]        size factor of SOI
%   tr_dotted[1]      [ 0 1 ] switch for a complete dotted transfer orbit
%
% OUTPUT:
%   r_SOI_1[1]        radius of SOI of the planet_1 [km]
%   r_SOI_2[1]        radius of SOI of the planet_2 [km] 
%   e_t[1]            eccentricity of the second transfer orbit [-] 
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
[kep_2_arr(:),~] = uplanet( T_t_2, p2 );
a_p2 = kep_2_arr(1);               % [km] semimajor axis Planet_1

Ratio_s1 = mean_radius_S / mr_p1;  % ratio between Planet_1's radius and Planet_2's radius
Ratio_12 = mr_p1 / mr_p2;          % ratio between Planet_1's radius and Planet_2's radius

r_SOI_1 = norm(R_1*AU)*(mu_p1/mu_S)^(2/5); % [km] radius of the Sphere of Influence
r_SOI_2 = norm(R_2*AU)*(mu_p2/mu_S)^(2/5); % [km] radius of the Sphere of Influence

%% Orbit
% Set options
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

if( nargin < 17 )
   tr_dotted = 0;end

% Set transfer trajectory data
Dt = T_t_2 - T_t_1;                % [days] duration of the transfer orbit
tspan_t = [ 0 Dt*3600*24 ];        % arc of ellipse of transfer orbit
[a_t,~,e_t,~,v_t_1,v_t_2,~,~] = lambertMR( r_t_1, r_t_2, Dt*24*3600 , mu_S, 0, 0, 0, 1);
if( e_t<1 && e_t>0 && tr_dotted==1 )
    T_t = (2*pi*sqrt(a_t^3/mu_S));     % [sec] period of the transfer orbit
    tspan_t_dot = [ 0 T_t-tspan_t(2)]; % dotted arc of ellipse
end
y0_t     = [ r_t_1; v_t_1' ];      % Initial conditions of the transfer orbit
y0_t_dot = [ r_t_2; v_t_2' ];      % Final conditions of the transfer orbit

% Set planets data
T_p1 = (2*pi*sqrt(a_p1^3/mu_S));  % [sec] period of Planet_1
T_p2 = (2*pi*sqrt(a_p2^3/mu_S));  % [sec] period of Planet_2
tspan_p1 = [ 0 T_p1 ];
tspan_p2 = [ 0 T_p2 ];
y0_p1 = [ r_t_1; v_p1 ];          % Initial conditions of Planet_1
y0_p2 = [ r_t_2; v_p2 ];          % Initial conditions of Planet_2

% Perturbation of the orbits
[~, Y_1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p1, y0_p1, options);
[~, Y_2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p2, y0_p2, options);
[~, Y_t] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t, y0_t, options);
if( e_t<1 && e_t>0 && tr_dotted==1 )
[~, Y_t_dot] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t_dot, y0_t_dot, options);end
%evaluate where is the Planet 1 at arrival time
[kep_1_arr(:),~] = uplanet( T_t_2, p1 );
    [ r_1_arr(:),~ ] = par2car( kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6), mu_S );
%evaluate where is the Planet 2 at departure time
[kep_2_arr(:),~] = uplanet( T_t_1, p2 );
    [ r_2_dep(:),~] = par2car( kep_2_arr(1), kep_2_arr(2), kep_2_arr(3), kep_2_arr(4), kep_2_arr(5), kep_2_arr(6), mu_S );

%% Plot Orbit Opt

if( nargin < 13 ) % with no size factors the plot will be correctly scaled
    Sun_fct = 1;
    p1_fct  = 1;
    p2_fct  = 1;
end
figure;hold on; grid on
% p_D = plot3 ( r_1(:,1), r_1(:,2), r_1(:,3), 'Color','#0072BD','LineWidth',15 ); %Departure window
% p_A = plot3 ( r_2(:,1), r_2(:,2), r_2(:,3), 'Color','#D95319','LineWidth',15 ); %Arrival window
% p_D.Color(4) = 0.3;
% p_A.Color(4) = 0.3;
plot3 ( Y_1(:,1),Y_1(:,2),Y_1(:,3), 'Color','#0072BD','LineWidth',2);
plot3 ( Y_2(:,1),Y_2(:,2),Y_2(:,3), 'Color','#D95319','LineWidth',2);
if( e_t<1 && e_t>0 && tr_dotted==1 )
plot3 ( Y_t_dot(:,1), Y_t_dot(:,2), Y_t_dot(:,3), 'Color','#77AC30','LineStyle','--','LineWidth',2);end
plot3 ( Y_t(:,1),Y_t(:,2),Y_t(:,3), 'Color','#77AC30','LineWidth',2);

factor_p1  = p1_fct * Sun_fct/Ratio_s1;                   % Size factor for planet_1
Planet_plot ( 10,[ 0 0 0 ], Sun_fct );                    % Sun plot
Planet_plot (  p1, r_t_1, factor_p1 );                    % Planet_1 plot departure
Planet_plot (  p1, r_1_arr, factor_p1 );                  % Planet_1 plot arrival
Planet_plot (  p2, r_t_2, p2_fct*factor_p1/Ratio_12 );    % Planet_2 plot arrival
Planet_plot (  p2, r_2_dep, p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot departure
% legend({'Departure window','Arrival window'});
% view(2);
view(35,20);
axis equal
if( nargin >= 16 ) % SOI plot
SOI_plot ( r_SOI_2*SOI_fct*Sun_fct, 0.5, r_t_2);end


s_sys_2.r_SOI_1 = r_SOI_1;
s_sys_2.r_SOI_2 = r_SOI_2;
s_sys_2.e_t = e_t;
end