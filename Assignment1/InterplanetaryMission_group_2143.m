% InterplanetaryMission_group_2143
%
% Main script for group 2143 for Assignment 1 - Interplanetary Mission
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
settings.timeWindows.plot = true;         % [-] True if plot are to be visulized
settings.timeWindows.parallel = true;     % [-] True if parallel computing is allowed

%%% OPTIMIZATION
data.optimization.minHfl = 200;

%%% GRID SEARCH
data.gridSearch.nDep  = 100;
data.gridSearch.nTOF1 = 100;
data.gridSearch.nTOF2 = 100;
settings.gridSearch.plot = true;
settings.gridSearch.parallel = true;

%% SELECTING TIME WINDOWS
data = timeWindows(data, settings);

%% GRID SEARCH OPTIMIZATION
data = gridSearch(data, settings);

%% PARTICLE SWARM OPTIMIZATION
[data] = optimizationPS(data, settings);

%% GA OPTIMIZATION
[data] = optimizationGA(data, settings);

%% FINAL OUTPUT
data.final.depDate = [2041, 8, 1, 3, 22, 17];
data.final.flybyDate = [2042, 6, 9, 10, 40, 51];
data.final.arrDate = [2042, 12, 6, 6, 59, 50];


plotMission(data);

















