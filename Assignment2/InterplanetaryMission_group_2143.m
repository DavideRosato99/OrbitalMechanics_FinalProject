%
% PlanetaryMission_group_2143
%
% Main script for group 2143 for Assignment 1
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
% -------------------------------------------------------------------------

close all
clear all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

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
%%% STARTING DATA
data.starting.depDate     = [2027 12 1 0 0 0];      % [Y M D H m S] Minimum departure date
data.starting.arrDate     = [2067 12 1 0 0 0];      % [Y M D H m S] Maximum arrival date
data.starting.depPlanID   = 4;                      % [-] Departure planet ID (Mars)
data.starting.flyByPlanID = 3;                      % [-] Fly-by planet ID (Earth)
data.starting.arrPlanID   = 2;                      % [-] Arrival planet ID (Venus)

%%% USED CONSTANTS
data.constants.muE = astroConstants(13);       % [km^3/s^2] Planetary constant of the Earth
data.constants.muS = astroConstants(4);        % [km^3/s^2] Planetary constant of the Sun
data.constants.AU  = astroConstants(2);        % [km] Astronomic unit
data.constants.Re  = astroConstants(23);       % [km] Earth radius

%%% TIME WINDOWS SELECTION
data.timeWindows.nDep  = 100;             % [days] delta days for departure date
data.timeWindows.nTOF1 = 800;             % [-] Number of points in which TOF1 will be calculated
data.timeWindows.nTOF2 = 100;             % [-] Number of points in which TOF2 will be calculated
settings.timeWindows.plot = false;         % [-] True if plot are to be visulized
settings.timeWindows.parallel = true;     % [-] True if parallel computing is allowed

%%% OPTIMIZATION
data.optimization.minHfl = 200;

%%% GRID SEARCH
data.gridSearch.nDep = 50;
data.gridSearch.nTOF1 = 50;
data.gridSearch.nTOF2 = 50;
settings.gridSearch.plot = true;
settings.gridSearch.parallel = true;

%% SELECTING TIME WINDOWS
% data = timeWindows(data, settings);

%% GRID SEARCH
data.timeWindows.depDate = [2039,7,8,7,42,58.9566981201172];
data.timeWindows.arrDate = [2055,7,29,22,53,46.4537894287109];
data.timeWindows.maxTOF1days = 2811.66781454150;
data.timeWindows.maxTOF2days = 3053.96467963780;
data = gridSearch(data, settings);


%% GRID SEARCH + FMINCON
% data.timeWindows.depDate = [2044,2,27,8,21,56.3843509521484];
% data.timeWindows.arrDate = [2055,6,18,14,24,49.3654221191406];
% [data] = gridSearch_fmincon(data, settings);

%% GA
% data.timeWindows.depDate = [2027 12 1 0 0 0];
% data.timeWindows.arrDate = [2067 12 1 0 0 0];
% [data] = optimizationGA(data, settings);























