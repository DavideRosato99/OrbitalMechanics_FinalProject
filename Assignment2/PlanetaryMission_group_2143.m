% PlanetaryMission_group_2143
%
% Main script for group 2143 for Assignment 2 - Planetary Mission
%
% CONTRIBUTORS:
%   Rosato Davide          10618468     davide.rosato@mail.polimi.it
%   Saba Mohammadi Yengeje 10789462     saba.mohammadi@mail.polimi.it
%   Spinelli Jason         10618465     jason.spinelli@mail.polimi.it
%   Tagliati Alessia       10635119     alessia.tagliati@mail.polimi.it
%
% NOTE: due to high computational cost of the codes, parallel computing is
%       preferred.
%
% VERSIONS
%   2021-10-21: Release
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
data.starting.OM = NaN;                         % [deg] Orbit RAAN
data.starting.om = NaN;                         % [deg] Orbit argument of pericenter
data.starting.th = 0;                          % [deg] Orbit true anomaly

%%% USED CONSTANTS
data.constants.muE = astroConstants(13);       % [km^3/s^2] Planetary constant of the Earth

%%% FIND OPTIMAL RAAN AND OM ..............................................
data.optimal.nPeriod = 13;                      % [-] Number of periods the orbit is propagated for
data.optimal.nPoints = 500;                    % [-] Number of points for each period at which coordinates are computed
data.optimal.nOM = 100;                         % [-] Number of RAAN that are calculated to find the optimal value
data.optimal.nom = 100;                         % [-] Number of omega that are calculated to find the optimal value
settings.optimal.plot = true;                  % [-] True if plot are to be visulized
settings.optimal.movie = true;                 % [-] True if movies are to be created
settings.optimal.parallel = true;              % [-] True if parallel computing is allowed

%%% GROUNDTRACKS
Tperiod = 2*pi * sqrt(data.starting.a^3/data.constants.muE);
data.groundtracks.periods = [Tperiod, 24*60*60, 10*24*60*60];   % [s] Periods for which the groundtracks will be displayed
data.groundtracks.k = 2;                                             % [-] Number of periods of the Earth
data.groundtracks.m = 5;                                             % [-] Number of periods of the satellite
settings.groundtracks.plot = true;                                   % [-] True if plot are to be visulized
settings.groundtracks.movie = true;                                  % [-] True if movies are to be created
settings.groundtracks.parallel = true;                               % [-] True if parallel computing is allowed

%%% PROPAGATION
data.propagation.TmaxComp = 40*Tperiod;               % [s] Final ode integration time for computational comparison
data.propagation.deltaSpan1 = 13000 : 1000 : 17000;   % [-] CPU settling number of N step vector
data.propagation.deltaSpan2 = 17000 : 5000 : 117000;  % [-] Computational comparison number of N step
data.propagation.Tmax = 10*365*24*3600;              % [s] Final ode integration time for orbit propagation
settings.propagation.plot = true;                     % [-] True if propagation plot are to be visualized
settings.propagation.movie = true;                    % [-] True if movies are to be created
settings.propagation.parallel = true;                 % [-] True if parallel computing is allowed

%%% REAL DATA COMPARISON
data.realData.Tmax = 2*365*24*60*60;
data.realData.nPoints = 365*200;


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
data = GroundTrack(data, settings);

%% PROPAGATE PERTURBED ORBIT WITH CARTESIAN AND GAUSS EQUATIONS
data = propagate(data, settings);

%% REAL DATA
load(strcat(pwd,'\functions\initialize\TLEs.mat'));
data = realData(data, satData, settings);

