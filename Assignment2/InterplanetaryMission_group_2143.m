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
data.timeWindows.nDep = 400;             % [days] delta days for departure date
data.timeWindows.nTOF1 = 400;              % [-] Number of points in which TOF1 will be calculated
data.timeWindows.nTOF2 = 3;              % [-] Number of points in which TOF2 will be calculated
settings.timeWindows.plot = true;         % [-] True if plot are to be visulized
settings.timeWindows.parallel = true;     % [-] True if parallel computing is allowed

%% SELECTING TIME WINDOWS
% [DVfirstLeg, DVsecondLeg, DVTOT, daysSpan, TOF1span, TOF2span, DVTOTreal, dtTOTreal]...
%     = timeWindows(depPlanID, flyByPlanID, arrPlanID, depDate,...
%     arrDate, deltaD);

data = prova(data, settings);



% %% PLOT
% figure
% I = 0 : 20 : (datenum(arrDate) - datenum(depDate));
% J = dtTOTreal;
% [a, b] = find(DVTOTreal ~= 0);
% [I, J] = meshgrid(I(a), J(b));
% 
% surf(DVTOTreal')
% 
% %%
% clc



%% FLY-BY

data.flyby.hAtm    = 200; %(?)
data.flyby.VVminus = %;
data.flyby.VVplus  = %;
 
[r_p, h_ga, Delta, delta_V_poweredFB,e_minus,e_plus,a_minus,a_plus] = ...
    our_powered_Flyby(data, settings)




