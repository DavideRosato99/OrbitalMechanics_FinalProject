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


%% SETTINGS
%%% FIND OPTIMAL RAAN AND OM
date0 = [2021 12 25 12 0 0];             % [-] Date time of the starting of orbit propagation
nPeriod = 10;                            % [-] Number of periods the orbit is propagated for
deltaPeriod = 1*24*60*60;                % [s] Delta period for TLEs parsing
nPoints = 100;                           % [-] Number of points for each period at which coordinates are computed
nOM = 36;                                % [-] Number of RAAN that are calculated to find the optimal value
nom = 36;                                % [-] Number of omega that are calculated to find the optimal value
plotOptimal = true;                      % [-] True if plot are to be visulized

%% USED CONSTANTS
muE = astroConstants(13);                % [km^3/s^2] Earth's gravitational parameter

%% INITIAL ORBIT
a = 7.2776e4;                 % [km] Orbit semi-major axis
e = 0.6665;                   % [-] Orbit eccentricity
i = 134.2783;                 % [deg] Orbit inclination
orbIn = [a, e, deg2rad(i)];

% Orbit period
Torbit = 2*pi * sqrt(a^3/muE);


%% RETRIEVE TLES
if exist(strcat(pwd, '\functions\initialize\TLEs.mat'), 'file')
    load(strcat(pwd, '\functions\initialize\TLEs.mat'));
    answer = questdlg({strcat("Actual date:              ", datestr(datetime('now'))),...
        strcat("Last TLEs update:    ", satData.lastUpdate), ...
        'Do you want to upload new TLEs or use the actual ones?'}, ...
        'Dialog', 'Existing', 'New', 'Abort', 'Abort');
    
    % handle response
    switch answer
        case 'Existing'
            fprintf('Using the TLEs file located in folder initialize... \n\n')
            load(strcat(pwd, '\functions\initialize\TLEs.mat'));
        case 'New'
            fprintf('Generating new TLEs...\n');
            retrieveTLEs
            fprintf('New TLEs file generated!\n\n');
            clearvars -except satData date0 nPeriod deltaPeriod nPoints nOM nom muE a e i orbIn Torbit
        case 'Abort'
            error('Simulation aborted')
    end
else
    fprintf('No TLEs file was found in folder initialize, generating a new one...\n\n');
    retrieveTLEs
end


%% FIN OPTIMAL RAAN AND OM ------------------------------------------------
% The following sections are performed in order to find the optimal OM and
% om which allows the satellite to follow a safe orbit, which points are as
% away as possible from other objects orbiting the Earth. Objects and
% satellite orbits are propagated for a total time of
% totTime = nPeriod*Tperiod, where:
%    - nPeriod: user-defined number of periods [-]
%    - Tperiod: given satellite orbit period   [s]
% Satellite orbit is propagated through an unperturbed model, while for the
% debris and other satellites orbiting the Earth NORAD TLEs are used. Since
% TLEs give reliable estimations for a maximum of 30 days from the day they
% are calculated, an error message is displayed if totTime is greater than
% 30 days. For computational efficiency purposes, only TLEs for those
% objects which have an orbit period which differs from unperturbed
% satellite orbit period of a maximum of deltaPeriod are computed.

%%% RETRIEVE TLEs DATA ----------------------------------------------------
% Get all earth orbiting objects positions and velocities (TEME reference frame)
satData = satData.data;
[rteme, vteme] = SGP4(satData, date0);
% Get TLEs objects orbital period to be compared to satellite's one
N = size(satData, 1);
T_TLEs = zeros(N, 1);

for i = 1:N
    [orb] = car2par(rteme(i, :), vteme(i, :), muE);
    T_TLEs(i) = 2*pi * sqrt(orb(1)^3/muE);
end

% Parsing TLEs
indexes = abs(T_TLEs - Torbit) < deltaPeriod;
raMax = 0;
ctrInd = 0;
for i = 1:N
    if indexes(i)
        ctrInd = ctrInd+1;
        [orb] = car2par(rteme(i, :), vteme(i, :), muE);
        ra = orb(1)*(1+orb(2));
        if ra > raMax
            raMax = ra;
        end
    end
end

%%% PROPAGATION OF TLEs OBJECTS -------------------------------------------
[Y, MO, D] = ymd(datetime(date0));
[H, M, S] = hms(datetime(date0));

% Define timespan
timespan = linspace(0, nPeriod*Torbit, nPeriod*nPoints);
RtleTEME = zeros(ctrInd, 3, length(timespan));
RtleECI = zeros(ctrInd, 3, length(timespan));
for i = 1:length(timespan)
    t = timespan(i);
    date = [Y MO D H M S+t];
    [Yl, MOl, Dl] = ymd(datetime(date));
    [Hl, Ml, Sl] = hms(datetime(date));
    % Get TLEs coordinates using SGP4 propagation algorithm
    [RtleTEME(:, :, i), ~] = SGP4(satData(indexes, :), [Yl MOl Dl Hl Ml Sl]);
    [RtleECI(:, :, i), ~] = teme2eci(RtleTEME(:, :, i), 0, date2mjd2000([Yl MOl Dl Hl Ml Sl]));
end


%%% SATELLITE UNPERTURBED ORBIT PROPAGATION -------------------------------
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
OmLoop  = deg2rad(linspace(1, 360, nOM));
omLoop  = deg2rad(linspace(1, 360, nom));
I = length(OmLoop);
J = length(omLoop);
Yorb = cell(I*J, 1);

% Propagato the orbit0
orb = [orbIn, 0, 0, 0];
[rr, vv] = par2car(orb, muE);
x0 = [rr; vv];
[T0, Y0] = ode113(@ode_2bp, timespan, x0, options, muE, 'cart');

R_i = [1 0 0; 0 cos(orbIn(3)) sin(orbIn(3)); 0 -sin(orbIn(3)) cos(orbIn(3))];
ctr = 1;
for i = 1:I
    Om = OmLoop(i);
    R_OM = [cos(Om) sin(Om) 0; -sin(Om) cos(Om) 0; 0 0 1];
    for j = 1:J
        om = omLoop(j);
        R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
        R = R_om * R_i * R_OM;
        Yorb{ctr} = (R*Y0(:, 1:3)')';
        ctr = ctr+1;
    end
end


%%% FIND MINIMUM DISTANCES ------------------------------------------------
MINdis = zeros(I, J);
MINdisMod = zeros(I, J);
ctr = 1;
for i = 1:I
    for j = 1:J
        minDis = zeros(ctrInd, 1);
        minDisMod = zeros(ctrInd, 1);
        for m = 1:ctrInd
            [minDis(m), ind] = min(vecnorm(Yorb{ctr}(:, 1:3)' - squeeze(RtleECI(m, :, :))));
            minDisMod(m) = minDis(m) * (1-exp(-T0(ind(1))/10000));
        end
        MINdis(i, j) = min(minDis);
        MINdisMod(i, j) = min(minDisMod);
        ctr = ctr + 1;
    end
end


%%% OPTIMAL ORBIT
BESTmod = max(max(MINdisMod))
[i1, i2] = find(MINdisMod == BESTmod)
OmLoop(i1)
omLoop(i2)
orb = [orbIn, OmLoop(i1), omLoop(i2), 0]

figure('Name','Satellite distance from TLEs - exponential','NumberTitle','off');
[OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
surf(OMM, OM, MINdisMod'); hold on;
shading interp
xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$f(\Omega,\omega)$')
BESTmod = max(max(MINdisMod));
[i1, i2] = find(MINdisMod == BESTmod);
plot3(OmLoop(i1)*180/pi, omLoop(i2)*180/pi, BESTmod, 'ro')
fprintf('\nOPTIMAL OM and om with exponential cost function:\nOM = %.2f\nom = %.2f\n\n', OmLoop(i1)*180/pi, omLoop(i2)*180/pi)


%% PROPAGATE UNPERTURBED ORBIT
% Time for: 1 orbit, 1 day, 10 days
Tvec = [Torbit, 24*60*60, 20*Torbit];
N = length(Tvec);
T = cell(N, 1);
Y = cell(N, 1);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];
labels = {'One orbit', 'One day', '10 days'};

figure
[Xs, Ys, Zs] = sphere(100);
Xs = 3671*Xs;
Ys = 3671*Ys;
Zs = 3671*Zs;
for i = 1:length(Tvec)
    [T{i}, Y{i}] = ode113(@ode_2bp, [0 Tvec(i)], Y0, options, muE, 'cart');
    subplot(1, 3, i)
    plot3(Y{i}(:,1), Y{i}(:,2), Y{i}(:,3)); hold on; axis equal; grid on
    surf(Xs, Ys, Zs)
    title(labels(i))
end

% %%% GROUNDTRACK
% GroundTrack(Tvec(3), orb, date0, 'unpert');
% 
% %%% REPEATING GROUNDTRACK
% k = 12;
% m = 1;
% GroundTrack(Tvec(3), orb, date0, 'unpert', m, k);

%% PROPAGATE PERTURBED ORBIT
% Time for: 1 orbit, 1 day, 10 days
Tvec = [Torbit, 24*60*60, 20*Torbit];
N = length(Tvec);
T = cell(N, 1);
Y = cell(N, 1);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];
labels = {'One orbit', 'One day', '10 days'};

figure
[Xs, Ys, Zs] = sphere(100);
Xs = 3671*Xs;
Ys = 3671*Ys;
Zs = 3671*Zs;

[Ye0, MO0, D0] = ymd(date0);
[H0, M0, S0] = hms(date0);

for i = 1:length(Tvec)
    [T{i}, Y{i}] = ode113(@ode_2bp, [0 Tvec(i)], Y0, options, muE, 'cart', ...
        datetime([Ye0, MO0, D0, H0, M0, S0]));
    subplot(1, 3, i)
    plot3(Y{i}(:,1), Y{i}(:,2), Y{i}(:,3)); hold on; axis equal; grid on
    surf(Xs, Ys, Zs)
    title(labels(i))
end

% %%% GROUNDTRACK
% GroundTrack(Tvec(3), orb, date0, 'pert');
% 
% %%% REPEATING GROUNDTRACK
% k = 12;
% m = 1;
% GroundTrack(Tvec(3), orb, date0, 'pert', m, k);



%% CARTESIAN AND GAUSS UNPERTURBED
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

%% CARTESIAN AND GAUSS PERTURBED
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];

[Ye0, MO0, D0] = ymd(date0);
[H0, M0, S0] = hms(date0);

[Tcart, Ycart] = ode113(@ode_2bp, [0 100*Torbit], Y0, options, muE, 'cart',...
    datetime([Ye0, MO0, D0, H0, M0, S0]));
[Tgauss, Ygauss] = ode113(@ode_2bp, [0 100*Torbit], orb, options, muE, 'gauss', ...
    datetime([Ye0, MO0, D0, H0, M0, S0]));

% for i = 1:length(Tgauss)
%     [Ygauss(i, 1:3), Ygauss(i, 4:6)] = par2car(Ygauss(i, :), muE);
% end

%%
for i = 1:length(Tcart)
    orb = car2par(Ycart(i,1:3), Ycart(i,4:6), muE);
    a(i) = orb(1);
end

figure
plot3(Ygauss(:,1),Ygauss(:,2),Ygauss(:,3))
legend('cart','gauss')


%% FILTER
[Tgauss, Ygauss] = ode113(@ode_2bp, [0 1000*Torbit], orb, options, muE, 'gauss', ...
    datetime([Ye0, MO0, D0, H0, M0, S0]));
[perturbations] = recallOdeFcn(@ode_2bp, Tgauss, Ygauss, muE, 'gauss', datetime([Ye0, MO0, D0, H0, M0, S0]));

rr = zeros(3, length(Tgauss)); vv = zeros(3, length(Tgauss));
for i = 1:length(Tgauss)
    [rr(:,i), vv(:,i)] = par2car(Ygauss(i,:), muE);
end

%%
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

%%
figure
colormap = parula(length(Tgauss));
for i = 2:10:length(Tgauss)
    plot3([rr(1,i) rr(1, i-1)], [rr(2,i) rr(2, i-1)], [rr(3,i) rr(3, i-1)], ...
        'color', colormap(i, :));
    hold on
end




%% PLOT ------------------------------------------------------------------
% figure('Name','Satellite distance from TLEs','NumberTitle','off');
% [OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
% surf(OMM, OM, MINdis'); hold on;
% shading interp
% xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$|r_{sat}-r_{TLEs}|_{min}$')
% BEST = max(max(MINdis));
% [i1, i2] = find(MINdis == BEST);
% plot3(OmLoop(i1)*180/pi, omLoop(i2)*180/pi, BEST, 'ro')
% fprintf('\nOPTIMAL OM and om:\nOM = %.2f\nom = %.2f\n\n', OmLoop(i1)*180/pi, omLoop(i2)*180/pi)
% 
% figure('Name','Satellite distance from TLEs - exponential','NumberTitle','off');
% [OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
% surf(OMM, OM, MINdisMod'); hold on;
% shading interp
% xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$f(\Omega,\omega)$')
% BESTmod = max(max(MINdisMod));
% [i1, i2] = find(MINdisMod == BESTmod);
% plot3(OmLoop(i1)*180/pi, omLoop(i2)*180/pi, BESTmod, 'ro')
% fprintf('\nOPTIMAL OM and om with exponential cost function:\nOM = %.2f\nom = %.2f\n\n', OmLoop(i1)*180/pi, omLoop(i2)*180/pi)
% 
% Yorbit = Yorb{(i1-1)*J + i2};
% 
% 
% %%
% figure('Name','Best orbit evolution','NumberTitle','off');
% [rteme, vteme] = SGP8(Y, D + hms2fracday(H, M, S), NORAD_TLEs(indexes,:));
% p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 1);
% axis equal; hold on; grid on
% plot3(Yorbit(:,1), Yorbit(:,2), Yorbit(:,3), 'r')
% sat = plot3(Yorbit(1,1), Yorbit(1,2), Yorbit(1,3), 'ro', 'MarkerSize', 3);
% xlim([-raMax raMax])
% ylim([-raMax raMax])
% zlim([-raMax raMax])
% [Xs, Ys, Zs] = sphere(100);
% Xs = 3671*Xs;
% Ys = 3671*Ys;
% Zs = 3671*Zs;
% surf(Xs, Ys, Zs)
% 
% for i = 1:length(timespan)
%     t = timespan(i);
%     date = datetime([Y MO D H M S+t]);
%     [Hl, Ml, Sl] = hms(date);
%     fracDay = day(date) + hms2fracday(Hl, Ml, Sl);
%     [rteme, vteme] = SGP8(year(date), fracDay, NORAD_TLEs(indexes,:));
%     delete(p)
%     delete(sat)
%     p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 2);
%     sat = plot3(Yorbit(i,1), Yorbit(i,2), Yorbit(i,3), 'ro', 'MarkerSize', 3);
%     drawnow limitrate
% end














