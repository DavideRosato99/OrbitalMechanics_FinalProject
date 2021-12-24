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


%% STARTING PARAMETERS
depDate = [2027 12 1 0 0 0];      % [-] Minimum departure date
arrDate = [2067 12 1 23 59 59];   % [-] Maximum arrival date
depPlanID   = 4;                  % [-] Departure planet ID
flyByPlanID = 3;                  % [-] Fly-by planet ID
arrPlanID   = 2;                  % [-] Arrival planet ID

%% SELECTING TIME WINDOWS
deltaD = 20;
[DVfirstLeg, DVsecondLeg, DVTOT, daysSpan, TOF1span, TOF2span, DVTOTreal, dtTOTreal]...
    = timeWindows(depPlanID, flyByPlanID, arrPlanID, depDate,...
    arrDate, deltaD);


%% PLOT
figure
I = 0 : 20 : (datenum(arrDate) - datenum(depDate));
J = dtTOTreal;
[a, b] = find(DVTOTreal ~= 0);
[I, J] = meshgrid(I(a), J(b));

surf(DVTOTreal')

%%
clc






