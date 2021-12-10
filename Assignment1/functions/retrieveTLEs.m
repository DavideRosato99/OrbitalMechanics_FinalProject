% RetrieveTLEs
%
% Script to retrieve and parse TLEs data from celestrak.com. It saves a
% .mat files containing several TLEs for different space objects. The
% file data.mat contains a table which sorts satellites TLEs through:
%     - name              char     Satellite name                    [-]
%     - epoch             double   Epoch                             [year]
%     - noradNumber       double   Catalog number (NORAD)            [-]
%     - bulletinNumber    double   Identification Number             [-]
%     - classification    char     Security Classification           [-]
%     - revolutionNumber  double   Revolution Number at Epoch        [-]
%     - ephemerisType     char     Ephemeris type                    [-]
%     - xm0               double   Mean Anomaly                      [rad]
%     - xnode0            double   Right Ascension of Ascending Node [rad]
%     - omega0            double   Argument of Perigee               [rad]
%     - xincl             double   Orbit Inlcination                 [rad]
%     - e0                double   Orbit Eccentricity                [-]
%     - xn0               double   Orbits per minute                 [1/min]
%     - xndt20            double   First time derivative of m.m.     [-]
%     - xndd60            double   Second time derivative of m.m.    [-]
%     - bstar             double   Bstar/drag term                   [-]
%
% NOTE: Seven different types of space satellites are retrieved and parsed:
%
%      1) Space Stations and active satellites
%      URL: https://www.celestrak.com/NORAD/elements/active.txt
%
%      2) Analyst Satellites 
%      URL: https://www.celestrak.com/NORAD/elements/analyst.txt
%
%      3) Russian ASAT Test Debris (COSMOS 1408)
%      URL: https://www.celestrak.com/NORAD/elements/1982-092.txt
%
%      4) Indian ASAT Test Debris
%      URL: https://www.celestrak.com/NORAD/elements/2019-006.txt
%
%      5) Fengyun 1C Debris
%      URL: https://www.celestrak.com/NORAD/elements/1999-025.txt
%
%      6) Iridium 33 Debris
%      URL: https://www.celestrak.com/NORAD/elements/iridium-33-debris.txt
%
%      7) Cosmos 2251 Debris
%      URL: https://www.celestrak.com/NORAD/elements/cosmos-2251-debris.txt
%
% ADVISE: The user is supposed to check the site 
%         https://www.celestrak.com/NORAD/elements/ and control if other
%         types of data are added.
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

clear all
close all
clc


%% SET THE URL FROM WHERE RETRIEVE
url = {'https://www.celestrak.com/NORAD/elements/active.txt'
       'https://www.celestrak.com/NORAD/elements/analyst.txt'
       'https://www.celestrak.com/NORAD/elements/1982-092.txt'
       'https://www.celestrak.com/NORAD/elements/2019-006.txt'
       'https://www.celestrak.com/NORAD/elements/1999-025.txt'
       'https://www.celestrak.com/NORAD/elements/iridium-33-debris.txt'
       'https://www.celestrak.com/NORAD/elements/cosmos-2251-debris.txt'};
   
   
%% INITIALIZATION
Name = '';
Cnum = zeros(1, 1);
SC = '';
ID = '';
year = zeros(1, 1);
doy = zeros(1, 1);
epoch = zeros(1, 1);
TD1 = zeros(1, 1);
TD2 = zeros(1, 1);
ExTD2 = zeros(1, 1);
BStar = zeros(1, 1);
ExBStar = zeros(1, 1);
Bstar = zeros(1, 1);
Etype = '';
Enum = zeros(1, 1);
incl = zeros(1, 1);
raan = zeros(1, 1);
e = zeros(1, 1);
omega = zeros(1, 1);
M = zeros(1, 1);
no = zeros(1, 1);
a = zeros(1, 1);
rNo = zeros(1, 1);


%% RETRIEVE DATA
N = size(url, 1);
for i = 1:N
    % Read and divide into lines data
    rawData = splitlines(webread(url{i}));
    M = size(rawData, 1);
    ctr = length(Name);
    for j = 1:M
        if not(isempty(strtrim(rawData{j})))
            
            if mod(j-1, 3) == 0 % Retrieve the name of the satellite
                ctr = ctr + 1;
                Name{ctr} = strtrim(rawData{j});
                
            elseif mod(j-2, 3) == 0 % Retrieve the first TLE line of the satellite
                line1 = strtrim(rawData{j});
                Cnum(ctr) = str2double(line1(3:7));
                SC{ctr} = line1(8);
                ID{ctr} = strtrim(line1(10:17));
                epoch(ctr) = str2double(line1(19:32));
                TD1(ctr) = str2double(line1(34:43));
                TD2(ctr) = str2double(line1(45:50));
                BStar(ctr) = str2double(line1(54:59));
                ExBStar(ctr) = str2double(line1(60:61));
                BStar(ctr) = BStar(ctr) * 1e-5 * 10^ExBStar(ctr);
                Etype{ctr} = strtrim(line1(63));
                
            elseif mod(j-3, 3) == 0 % Retrieve the second TLE line of the satellite
                line2 = strtrim(rawData{j});
                incl(ctr) = str2double(line2(9:16));
                raan(ctr) = str2double(line2(18:25));
                e(ctr) = str2double(strcat('0.', line2(27:33)));
                omega(ctr) = str2double(line2(35:42));
                M(ctr) = str2double(line2(44:51));
                no(ctr) = str2double(line2(53:63));
                rNo(ctr) = str2double(line2(65:end));
                
            end
        end
    end
end


%% SAVE DATA
NORAD_TLEs = table(Name', epoch', Cnum', ID', SC', rNo', Etype', M'*(pi/180), ...
    raan'*(pi/180), omega'*(pi/180), incl'*(pi/180), e', no'*(2*pi/1440), ...
    TD1' * 1e-8 * (2*pi/(1440^2)), TD2' * (2*pi/(1440^3)), BStar');
NORAD_TLEs.Properties.VariableNames = {'name', 'epoch', 'noradNumber', ...
    'bulletinNumber', 'classification', 'revolutionNumber', ...
    'ephemerisType', 'xm0', 'xnode0', 'omega0', 'xincl', 'e0', 'xn0', ...
    'xndt20', 'xndd60', 'bstar'};

save(strcat(pwd, '/functions/initialize/NORAD_TLEs'), 'NORAD_TLEs')
clearvars



