% PlanetaryMission_group_2143
%
% Main script for group 2143 for Assignment 1 - Planetary Explorer Mission
%
% CONTRIBUTORS:
%   Rosato Davide          10618468     davide.rosato@mail.polimi.it
%   Saba Mohammadi Yengeje 10789462     saba.mohammadi@mail.polimi.it
%   Spinelli jason         10618465     jason.spinelli@mail.polimi.it
%   Tagliati Alessia       10635119     alessia.tagliati@mail.polimi.it
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

close all
clear all
clc
set(0, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultTextInterpreter', 'latex')


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


%% GEt TLEs
% the following are performed in order to retrieve last-updated TLEs of
% satellites, space stations and debris orbiting around Earth
if not(exist(strcat(pwd, '/functions/initialize/NORAD_TLEs.mat'), 'file'))
    retrieveTLEs;
end
load(strcat(pwd, '/functions/initialize/NORAD_TLEs.mat'))


%% SETTINGS
%%% FIND OPTIMAL RAAN AND OM
date = datetime([2021 12 25 12 0 0]);       % [-] Date time of the starting of orbit propagation
nPeriod = 10;                 % [-] Number of periods the orbit is propagated for
deltaPeriod = 1*24*60*60;       % [s] Delta period for TLEs parsing
nPoints = 100;                % [-] Number of points for each period at which coordinates are computed
nOM = 30;                     % [-] Number of RAAN that are calculated to find the optimal value
nom = 30;                     % [-] Number of omega that are calculated to find the optimal value


%% USED CONSTANTS
muE = astroConstants(13);     % [km^3/s^2] Earth's gravitational parameter


%% INITIAL ORBIT
a = 7.2776e4;                 % [km] Orbit semi-major axis
e = 0.6665;                   % [-] Orbit eccentricity
i = 134.2783;                 % [deg] Orbit inclination
orbIn = [a, e, deg2rad(i)];

% Orbit period
Torbit = 2*pi * sqrt(a^3/muE);


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
[rteme, vteme] = SGP8(0, 0, NORAD_TLEs);
% Get TLEs objects orbital period to be compared to satellite's one
N = size(NORAD_TLEs, 1);
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
Y = year(date);
MO = month(date);
D = day(date);
[H, M, S] = hms(date);

% Define timespan
timespan = linspace(0, nPeriod*Torbit, nPeriod*nPoints);
Rtle = zeros(ctrInd, 3, length(timespan));
for i = 1:length(timespan)
    t = timespan(i);
    date0 = [Y MO D H M S+t];
    Yearl = year(datetime(date0));
    Dayl = day(datetime(date0));
    [Hl, Ml, Sl] = hms(datetime(date0));
    % Get TLEs coordinates using SGP8 propagation algorithm
    [Rtle(:, :, i), ~] = SGP8(Yearl, Dayl + hms2fracday(Hl, Ml, Sl), NORAD_TLEs(indexes,:));
end


%%% SATELLITE UNPERTURBED ORBIT PROPAGATION -------------------------------
options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
OmLoop = deg2rad(linspace(1, 360, nOM));
omLoop = deg2rad(linspace(1, 360, nom));
I = length(OmLoop);
J = length(omLoop);
Yorb = cell(I*J, 1);

% Propagato the orbit0
orb = [orbIn, 0, 0, 0];
[rr, vv] = par2car(orb, muE);
x0 = [rr; vv];
[T0, Y0] = ode113(@ode_2bp, timespan, x0, options, muE);

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
            [minDis(m), ind] = min(vecnorm(Yorb{ctr}(:, 1:3)' - squeeze(Rtle(m, :, :))));
            minDisMod(m) = minDis(m) * (1-exp(-T0(ind(1))/10000));
        end
        MINdis(i, j) = min(minDis);
        MINdisMod(i, j) = min(minDisMod);
        ctr = ctr + 1;
    end
end


%% PLOT ------------------------------------------------------------------
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














