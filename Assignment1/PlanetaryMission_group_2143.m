%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             ASSIGNMENT 1                                %
%                       Planetary Explorer Mission                        %
%                                                                         %
%                              GROUP 2143                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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

%% CALCULATIONS
tsince = 0;
t = tic;
[rTLEs_TEME, vTLEs_TEME] = SGP8(tsince, NORAD_TLEs(1:500,:));
toc(t)

error('ciao')
figure
p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 1);
axis equal; hold on; grid on
[Xs, Ys, Zs] = sphere(100);
Xs = 3671*Xs;
Ys = 3671*Ys;
Zs = 3671*Zs;
surf(Xs, Ys, Zs)

for t = 0:0.1:1440
[rteme, vteme] = SGP8(t, NORAD_TLEs(1:500,:));
delete(p)
p = plot3(rteme(:, 1), rteme(:, 2), rteme(:, 3), 'go', 'MarkerSize', 1);
drawnow limitrate
end














