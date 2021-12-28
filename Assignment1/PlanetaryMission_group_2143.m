% PlanetaryMission_group_2143
%
% Main script for group 2143 for Assignment 1 - Planetary Explorer Mission
%
% CONTRIBUTORS:
%   Rosato Davide          10618468     davide.rosato@mail.polimi.it
%   Saba Mohammadi Yengeje 10789462     saba.mohammadi@mail.polimi.it
%   Spinelli Jason         10618465     jason.spinelli@mail.polimi.it
%   Tagliati Alessia       10635119     alessia.tagliati@mail.polimi.it
%
% VERSIONS
%   2021-10-21: Release
%
% Copyright © 2021, All rights reserved
% SPDX-License-Identifier: GPL-3.0-or-later
%
% -------------------------------------------------------------------------

close all
clear all
clc
set(0, 'defaultTextInterpreter', 'latex')


%% SET PATH
% global path
filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;
if not(strcmp(filePath, currentPath))
    cd (filePath);
    currentPath = filePath;
end

addpath(genpath(currentPath));

% functions path
addpath(genpath('functions'))

%% **************************** USER HERE *********************************
% ! ATTENTION ! if OM and om are set to NaN, the first part of the script,
% the choice of the latters, will be performed. Choose a finite value for
% OM and om to skip the first part.

%%% INITIAL ORBIT PARAMETERS ..............................................
data.starting.date = [2021 12 20 12 0 0];      % [-]   Date time of the starting of orbit propagation
data.starting.a = 7.2776e4;                    % [km]  Orbit semi-major axis
data.starting.e = 0.6665;                      % [-]   Orbit eccentricity
data.starting.i = 134.2783;                    % [deg] Orbit inclinationdata.starting.a = 7.2776e4;                    
data.starting.OM = 72;                         % [deg] Orbit RAAN
data.starting.om = 45;                         % [deg] Orbit argument of pericenter
data.starting.th = 0;                          % [deg] Orbit true anomaly

%%% USED CONSTANTS
data.constants.muE = astroConstants(13);       % [km^3/s^2] Planetary constant of the Earth

%%% FIND OPTIMAL RAAN AND OM ..............................................
data.optimal.nPeriod = 1;                      % [-] Number of periods the orbit is propagated for
data.optimal.nPoints = 100;                    % [-] Number of points for each period at which coordinates are computed
data.optimal.nOM = 50;                         % [-] Number of RAAN that are calculated to find the optimal value
data.optimal.nom = 50;                         % [-] Number of omega that are calculated to find the optimal value
settings.optimal.plot = true;                  % [-] True if plot are to be visulized
settings.optimal.movie = true;                 % [-] True if movies are to be created
settings.optimal.parallel = true;              % [-] True if parallel computing is allowed

%%% GROUNDTRACKS
Tperiod = 2*pi * sqrt(data.starting.a^3/data.constants.muE);
data.groundtracks.periods = [3.25*Tperiod, 24*60*60, 10*24*60*60];   % [s] Periods for which the groundtracks will be displayed
data.groundtracks.k = 2;                                             % [-] Number of periods of the Earth
data.groundtracks.m = 5;                                             % [-] Number of periods of the satellite
settings.groundtracks.plot = true;                                   % [-] True if plot are to be visulized
settings.groundtracks.movie = true;                                  % [-] True if movies are to be created
settings.groundtracks.parallel = true;                               % [-] True if parallel computing is allowed

%%% PROPAGATION
%data.propagation.TmaxComp = 40*Tperiod;               % [s] Final ode integration time for computational comparison
%data.propagation.deltaSpan1 = 13000 : 1000 : 17000;    % [-] CPU settling number of N step vector
%data.propagation.deltaSpan2 = 17000 : 5000 : 117000;   % [-] Computational comparison number of N step
data.propagation.TmaxComp = Tperiod;               % [s] Final ode integration time for computational comparison
data.propagation.deltaSpan1 = 13000;    % [-] CPU settling number of N step vector
data.propagation.deltaSpan2 = 17000;   % [-] Computational comparison number of N step
data.propagation.Tmax = Tperiod;              % [s] Final ode integration time for orbit propagation
settings.propagation.plot = true;                     % [-] True if propagation plot are to be visualized
settings.propagation.movie = true;                    % [-] True if movies are to be created
settings.propagation.parallel = true;                 % [-] True if parallel computing is allowed
%%% REAL DATA COMPARISON
data.realData.Tmax = 365*24*60*60;
data.realData.nPoints = 365*100;


%% **** FROM NOW ON, DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING *****
% Retrieve given data
a   = data.starting.a;
e   = data.starting.e;
i   = data.starting.i;
OM  = data.starting.OM;
om  = data.starting.om;
th  = data.starting.th;
muE = data.constants.muE;

% Calculation
data.starting.orbIn = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];        % [-] Given orbit initial parameters
data.starting.Torbit = 2*pi * sqrt(a^3/muE);                          % [s] Given orbit period

clearvars -except data settings

%% RETRIEVE TLES ----------------------------------------------------------
if isnan(data.starting.OM) && isnan(data.starting.om)
    if exist(strcat(pwd, '\functions\initialize\TLEs.mat'), 'file')
        load(strcat(pwd, '\functions\initialize\TLEs.mat'));
        answer = questdlg({"WARNING: Since TLEs are reliable for a maximum of 30 days, use a TLEs data file uploaded in a date which is as close as possible to the satellite departure date.",...
            "  ", strcat("Selected departure date:     ", ...
            datestr(datetime(data.starting.date))),...
            strcat("Last TLEs update:               ", satData.lastUpdate), ...
            strcat("Actual date:                         ", datestr(datetime('now'))),...
            "   ",...
            'Do you want to upload actual date TLEs from Internet or use the existing one?'}, ...
            'Dialog', 'Existing', 'New', 'Abort', 'Abort');

        % handle response
        switch answer
            case 'Existing'
                fprintf('Using the TLEs file located in folder initialize... \n\n')
                load(strcat(pwd, '\functions\initialize\TLEs.mat'));
            case 'New'
                fprintf('Generating new TLEs from Internet...\n');
                retrieveTLEs
                fprintf('New TLEs file generated!\n\n');
                clearvars -except data satData settings
            case 'Abort'
                error('Simulation aborted')
        end
    else
        fprintf('No TLEs file was found in folder initialize, generating a new one...\n\n');
        retrieveTLEs
    end
end


%% FIN OPTIMAL RAAN AND OM ------------------------------------------------
if isnan(data.starting.OM) && isnan(data.starting.om)
    data = optimalParam(data, satData, settings);
end


%% GROUNDTRACKS
%% data = GroundTrack(data, settings);


%% PROPAGATE PERTURBED ORBIT WITH CARTESIAN AND GAUSS EQUATIONS
%% data = propagate(data, settings);

%% REAL DATA
load(strcat(pwd,'\functions\initialize\TLEs.mat'));
data = realData(data, satData, settings);






%% PROPAGATE PERTURBED ORBIT
% %Time for: 1 orbit, 1 day, 10 days
% Tvec = [Torbit, 24*60*60, 20*Torbit];
% N = length(Tvec);
% T = cell(N, 1);
% Y = cell(N, 1);
% options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
% 
% [rr0, vv0] = par2car(orb, muE);
% Y0 = [rr0; vv0];
% labels = {'One orbit', 'One day', '10 days'};
% 
% figure
% [Xs, Ys, Zs] = sphere(100);
% Xs = 3671*Xs;
% Ys = 3671*Ys;
% Zs = 3671*Zs;
% 
% [Ye0, MO0, D0] = ymd(date0);
% [H0, M0, S0] = hms(date0);
% 
% for i = 1:length(Tvec)
%     [T{i}, Y{i}] = ode113(@ode_2bp, [0 Tvec(i)], Y0, options, muE, 'cart', ...
%         datetime([Ye0, MO0, D0, H0, M0, S0]));
%     subplot(1, 3, i)
%     plot3(Y{i}(:,1), Y{i}(:,2), Y{i}(:,3)); hold on; axis equal; grid on
%     surf(Xs, Ys, Zs)
%     title(labels(i))
% end

%%% GROUNDTRACK
GroundTrack(Tvec(3), orb, date0, 'pert');

%%% REPEATING GROUNDTRACK
k = 12;
m = 1;
GroundTrack(Tvec(3), orb, date0, 'pert', m, k);
% 
% 
% 
% CARTESIAN AND GAUSS UNPERTURBED
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];
[Tcart, Ycart] = ode113(@ode_2bp, [0 100*Torbit], Y0, options, muE, 'cart');
orb
[Tgauss, Ygauss] = ode113(@ode_2bp, [0 100*Torbit], orb, options, muE, 'gauss');

for i = 1:length(Tgauss)
    [Ygauss(i, 1:3), Ygauss(i, 4:6)] = par2car(Ygauss(i, :), muE);
end
figure
plot3(Ycart(:,1),Ycart(:,2),Ycart(:,3))
hold on
plot3(Ygauss(:,1),Ygauss(:,2),Ygauss(:,3))

% CARTESIAN AND GAUSS PERTURBED
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];

[Ye0, MO0, D0] = ymd(date0);
[H0, M0, S0] = hms(date0);

[Tcart, Ycart] = ode113(@ode_2bp, [0 100*Torbit], Y0, options, muE, 'cart',...
    datetime([Ye0, MO0, D0, H0, M0, S0]));
[Tgauss, Ygauss] = ode113(@ode_2bp, [0 100*Torbit], orb, options, muE, 'gauss', ...
    datetime([Ye0, MO0, D0, H0, M0, S0]));

for i = 1:length(Tgauss)
    [Ygauss(i, 1:3), Ygauss(i, 4:6)] = par2car(Ygauss(i, :), muE);
end

%
for i = 1:length(Tcart)
    orb = car2par(Ycart(i,1:3), Ycart(i,4:6), muE);
    a(i) = orb(1);
end

figure
plot3(Ygauss(:,1),Ygauss(:,2),Ygauss(:,3))
legend('cart','gauss')


% FILTER
[Tgauss, Ygauss] = ode113(@ode_2bp, [0 1000*Torbit], orb, options, muE, 'gauss', ...
    datetime([Ye0, MO0, D0, H0, M0, S0]));
[perturbations] = recallOdeFcn(@ode_2bp, Tgauss, Ygauss, muE, 'gauss', datetime([Ye0, MO0, D0, H0, M0, S0]));

rr = zeros(3, length(Tgauss)); vv = zeros(3, length(Tgauss));
for i = 1:length(Tgauss)
    [rr(:,i), vv(:,i)] = par2car(Ygauss(i,:), muE);
end

%
close all

figure('Name','a','NumberTitle','off');
plot(Tgauss, Ygauss(:,1)); grid on;
hold on;
plot(Tgauss, movmean(Ygauss(:,1), 15000), 'LineWidth', 2)
title('Semi major axis');
xlabel('Time [s]'); ylabel('a [km]');

figure('Name','e','NumberTitle','off');
plot(Tgauss, Ygauss(:,2)); grid on;
hold on;
plot(Tgauss, movmean(Ygauss(:,2), 15000), 'LineWidth', 2)
title('Eccentricity');
xlabel('Time [s]'); ylabel('e [-]');

figure('Name','i','NumberTitle','off');
plot(Tgauss, wrapTo360(Ygauss(:,3)*180/pi)); grid on;
hold on;
plot(Tgauss, wrapTo360(movmean(Ygauss(:,3)*180/pi, 15000)), 'LineWidth', 2)
title('Inclination');
xlabel('Time [s]'); ylabel('i [deg]');

figure('Name','aJ2','NumberTitle','off');
plot(Tgauss, vecnorm(perturbations.aJ2)); grid on;
hold on;
title('aJ2');
xlabel('Time [s]'); ylabel('aJ2 [km/$s^2$]');

figure('Name','aMOON','NumberTitle','off');
plot(Tgauss, vecnorm(perturbations.aMOON)); grid on;
hold on;
title('aMOON');
xlabel('Time [s]'); ylabel('aMOON [km/$s^2$]');

figure('Name','perturbations acc.','NumberTitle','off');
plot(Tgauss, vecnorm(perturbations.aMOON + perturbations.aJ2)); grid on;
hold on;
plot(Tgauss, movmean(vecnorm(perturbations.aMOON + perturbations.aJ2), Tgauss(end), 'SamplePoints', Tgauss,...
    'endpoints', 'shrink'), 'LineWidth', 2)
title('a Pert');
xlabel('Time [s]'); ylabel('a pert [km/$s^2$]');

figure('Name','Position','NumberTitle','off');
plot(Tgauss, vecnorm(rr)); grid on;
hold on;
plot(Tgauss, movmean(vecnorm(rr), Tgauss(end), 'SamplePoints', Tgauss,...
    'endpoints', 'shrink'), 'LineWidth', 2)
title('rr');
xlabel('Time [s]'); ylabel('r [km]');

%
figure
colormap = parula(length(Tgauss));
for i = 2:10:length(Tgauss)
    plot3([rr(1,i) rr(1, i-1)], [rr(2,i) rr(2, i-1)], [rr(3,i) rr(3, i-1)], ...
        'color', colormap(i, :));
    hold on
end




% PLOT ------------------------------------------------------------------
figure('Name','Satellite distance from TLEs','NumberTitle','off');
[OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
surf(OMM, OM, MINdis'); hold on;
shading interp
xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$|r_{sat}-r_{TLEs}|_{min}$')
BEST = max(max(MINdis));
[i1, i2] = find(MINdis == BEST);
plot3(OmLoop(i1)*180/pi, omLoop(i2)*180/pi, BEST, 'ro')
fprintf('\nOPTIMAL OM and om:\nOM = %.2f\nom = %.2f\n\n', OmLoop(i1)*180/pi, omLoop(i2)*180/pi)

figure('Name','Satellite distance from TLEs - exponential','NumberTitle','off');
[OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
surf(OMM, OM, MINdisMod'); hold on;
shading interp
xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$f(\Omega,\omega)$')
BESTmod = max(max(MINdisMod));
[i1, i2] = find(MINdisMod == BESTmod);
plot3(OmLoop(i1)*180/pi, omLoop(i2)*180/pi, BESTmod, 'ro')
fprintf('\nOPTIMAL OM and om with exponential cost function:\nOM = %.2f\nom = %.2f\n\n', OmLoop(i1)*180/pi, omLoop(i2)*180/pi)

Yorbit = Yorb{(i1-1)*J + i2};


%%
figure('Name','Best orbit evolution','NumberTitle','off');
[rteme, vteme] = SGP8(Y, D + hms2fracday(H, M, S), NORAD_TLEs(indexes,:));
p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 1);
axis equal; hold on; grid on
plot3(Yorbit(:,1), Yorbit(:,2), Yorbit(:,3), 'r')
sat = plot3(Yorbit(1,1), Yorbit(1,2), Yorbit(1,3), 'ro', 'MarkerSize', 3);
xlim([-raMax raMax])
ylim([-raMax raMax])
zlim([-raMax raMax])
[Xs, Ys, Zs] = sphere(100);
Xs = 3671*Xs;
Ys = 3671*Ys;
Zs = 3671*Zs;
surf(Xs, Ys, Zs)

for i = 1:length(timespan)
    t = timespan(i);
    date = datetime([Y MO D H M S+t]);
    [Hl, Ml, Sl] = hms(date);
    fracDay = day(date) + hms2fracday(Hl, Ml, Sl);
    [rteme, vteme] = SGP8(year(date), fracDay, NORAD_TLEs(indexes,:));
    delete(p)
    delete(sat)
    p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 2);
    sat = plot3(Yorbit(i,1), Yorbit(i,2), Yorbit(i,3), 'ro', 'MarkerSize', 3);
    drawnow limitrate
end












