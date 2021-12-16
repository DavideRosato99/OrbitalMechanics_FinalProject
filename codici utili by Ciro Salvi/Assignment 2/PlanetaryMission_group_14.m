close all; clear all; clc
% The script Explorer defines all necessary data for the assignment about
% the planetary explorer mission. It is divided into different sections:
%
% - Orbit Characterisics section evaluates the ground tracks for the
%   different conditions stated in the assignment and finds repeating GTs
% - Cartesian propagates the perturbed orbit with cartesian coordinates
% - Gauss propagates the orbit through Gauss equations
% - Comparison plots the keplerian elements found with Gauss and Cartesian,
%   for a selectable interval of time, and plots the error between the two
%   keplerian elements sets
% - Molniya 1-36 takes TLE data from similar-eccentricity and similar-semi
%   major axis satellite, launched in the late 70s
%
% All the sections can be switched on or off by setting to 1 or 0 the
% correspondant "check_..." switch. In all the sections it's possible to
% turn on/off different switches, later on specified
% 
% 
% CONTRIBUTORS
%   Ciro Salvi
%   Fabio Spaada
%   Alessandro Staffolani
%   Suhailah Alkhawashke
% 
% 

addpath('functions')
profile on
name = 'Planetary Explorer'; % name of the satellite
% --- Keplerian element ---
         
a  = 23210;                % [km] Semimajor axis
e  = 0.7091;               % [-] Eccentricity
i  = deg2rad(12.6865);     % [rad] Inclination
OM = deg2rad(0);           % [rad] Longitude of ascending node     4.9752;
om = deg2rad(0);           % [rad] Argument of periapsis           4.7063; 
f0 = deg2rad(0);           % [rad] True anomaly
k  = 12;                   % number of satellite's periods
m  = 5;                    % number of Earth's revolutions
hp = 6751.789;             % [km] 
cd = 2.1;                  % Drag coefficient
am_ratio = 0.060;          % [m^2/kg]

mu  = astroConstants(13);  % [km^3/s^2] Earth's planetary constants
R_e = astroConstants(23);  % [km] Earth's mean radius
J_2 = astroConstants(9); %0.00108263;          % second zonal harmonic 9 for astroConstants.m

w_E = 15.04;               % [deg/h] Earth' angular velocity
fg0 = deg2rad(30);

T_sat = 2*pi*sqrt(a^3/mu); % [s] Satellite's period
T_earth = 23*3600+56*60;   % sidereal Earth's time
time.T_sat   = duration(seconds(T_sat),'Format','dd:hh:mm:ss');
time.T_earth = duration(seconds(T_earth),'Format','dd:hh:mm:ss');


%% Orbit Characteristics
check_orbit_char = 0;
switch_plot = 'repeating_j2';  % 'all' || 'original' || 'repeating' || 'perturbed' || 'repeating_j2' || 'original_repeating' || 'repeating_perturbed' || 'original_perturbed'
                           % this switch plots the groud track in different 
                           % -above specified- conditions
time_evaluation = '12-revolution'; % 'one-day' || 'ten-days' || 'one-revolution' || '12-revolution'
graphs_a_initial_cond = 0;  % to plot the graphs for the initial condition search


if check_orbit_char
    check_orbit = 0;           % this switch plots the 3d orbit
if( strcmp(time_evaluation, 'one-day') )
    tspan = linspace(0,T_earth,10000);
    elseif( strcmp(time_evaluation, 'ten-days'))
    tspan = linspace(0,10*T_earth,10000);
    else
    tspan = NaN; % tspan will be evaluated fot each ground_track
end
    Orbit_Characteristics
end
if graphs_a_initial_cond && check_orbit_char
    try_GT_kep
end

moln_comparison = 0;        % compares the gauss propagation mathod with reality.

%% Cartesian
check_cartesian = 1;
if check_cartesian
J2_pert = 'on';             % switch 'on' to consider the J2 perturbation
drag_pert = 'on';           % switch 'on' to consider the drag perturbation
check_energy = 1;           % '1' to evaluate h,eps,e
check_graphs = 1;           % '1' to have their graphs
check_kep = 1;              % '1' to compute keplerian elements
check_graphs_kep = 1;       % '1' to haave keplerian elements graphs
check_draw = 1;             % '1' to plot the orbit
n_revolutions  = 12;        % number of periods related to the initial semi
                            % major axis considered.  
                            % if the satellite crashes, the propagation
                            % ends automatically when this condition occurs
                            
                            
if check_graphs
    check_energy = 1;
end
if check_graphs_kep
    check_kep = 1;
end
    Cartesian_Perturbed
end

%% Gauss
check_gauss = 1;

if check_gauss
J2_pert = 'on';             % switch 'on' to consider the J2 perturbation
drag_pert = 'on';           % switch 'on' to consider the drag perturbation
check_graphs_kep = 1;       % '1' to have keplerian elements graphs
check_energy = 1;           % '1' to evaluate h,eps,e
check_graphs = 1;           % '1' to have their graphs
check_draw = 1;             % '1' to plot the orbit
number_revolutions  = 1;   % same as above (Cartesian)



if check_graphs
    check_energy = 1;
end
    Gauss_perturbed
end

%% comparison
check_kep_comparison = 1;   
                            
        if check_cartesian == 0 || check_gauss == 0 
            check_kep_comparison = 0;
        end

if check_kep_comparison
time_start = 0;             % time from which the comparison is performed
number_of_periods = 5000;   % number of orbital periods to plot
check_plots = 1;            % '1' to plot singularly each kep with the error
check_subplots = 1;         % '1' to plot all the kep in the same figure

    comparison
end

%% molniya 1-36
check_molniya = 0;
if check_molniya
check_orbit = 1;            % '1' to plot the orbit wrt time
check_graphs_moln = 1;      % '1' to have keplerian elements evolution wrt time
moln_comparison = 0;        % compares the gauss propagation method with reality.



J2_pert = 'on';             % switch 'on' to consider the J2 perturbation
drag_pert = 'on';           % switch 'on' to consider the drag perturbation
check_graphs_kep = 1;       % '1' to have keplerian elements graphs
check_energy = 0;           % '1' to evaluate h,eps,e
check_graphs = 0;           % '1' to have their graphs
check_draw = 0;             % '1' to plot the orbit
number_revolutions  = 2455;   % same as above (Cartesian)

% time_start = 5750;
% Tcut = 6600;

% time_start = 2100;
% Tcut = 3300;
time_start = 0;
Tcut = 4000;

    ephemerides
if moln_comparison
    
    Gauss_perturbed
    figure()
        plot(T_Gauss/86400,filter_mov(kep_Gauss(:,1),T_Gauss, kep_Gauss,mu,3),Tmol_cut-timec_start,kepmol_cut(:,1))
        grid on
        xlabel ('time [days]')
        ylabel ('a [km]')
        title ('filtered semimajor axis evolution (and filtered)')
        legend('propagated and filtered','real')
    figure()
        plot(T_Gauss/86400,kep_Gauss(:,2),T_Gauss/86400,filter_mov(kep_Gauss(:,2),T_Gauss, kep_Gauss,mu,2),Tmol_cut-time_start,kepmol_cut(:,2))
        grid on
        xlabel ('time [days]')
        ylabel ('e')
        title ('eccentricity evolution (and filtered)')
        legend('propagated','filtered','real')
    figure() 
        plot(T_Gauss/86400,kep_Gauss(:,3)*180/pi,T_Gauss/86400,filter_mov(kep_Gauss(:,3)*180/pi,T_Gauss, kep_Gauss,mu,2),Tmol_cut-time_start,kepmol_cut(:,3)*180/pi)
        grid on
        xlabel ('time [days]')
        ylabel ('i [deg]')
        title ('inclination evolution (and filtered)')
        legend('propagated','filtered','real')
    figure()
        plot(T_Gauss/86400,kep_Gauss(:,4)*180/pi,Tmol_cut-time_start,kepmol_cut(:,4)*180/pi)
        grid on
        xlabel ('time [days]')
        ylabel ('OM [deg]')
        title ('RAAN evolution (and filtered)')
        legend('propagated','real')
    figure()
        plot(T_Gauss/86400,kep_Gauss(:,5)*180/pi,Tmol_cut-time_start,kepmol_cut(:,5)*180/pi)
        grid on
        xlabel ('time [days]')
        ylabel ('om [deg]')
        title ('argument of periapsis evolution (and filtered)')
        legend('propagated','real')
       

end
end

ppp = profile('info');
