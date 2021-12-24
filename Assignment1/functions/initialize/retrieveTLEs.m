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
%   Spinelli Jason         10618465     jason.spinelli@mail.polimi.it
%   Tagliati Alessia       10635119     alessia.tagliati@mail.polimi.it
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------


%% SET THE URL FROM WHERE RETRIEVE
url = {'https://www.celestrak.com/NORAD/elements/active.txt'
       'https://www.celestrak.com/NORAD/elements/analyst.txt'
       'https://www.celestrak.com/NORAD/elements/1982-092.txt'
       'https://www.celestrak.com/NORAD/elements/2019-006.txt'
       'https://www.celestrak.com/NORAD/elements/1999-025.txt'
       'https://www.celestrak.com/NORAD/elements/iridium-33-debris.txt'
       'https://www.celestrak.com/NORAD/elements/cosmos-2251-debris.txt'};
   
   
%% INITIALIZATION
satData = struct;
N = size(url, 1);
ctr = 0;

year = zeros(1, 1);
month = zeros(1, 1);
day = zeros(1, 1);
hr = zeros(1, 1);
minute = zeros(1, 1);
sec = zeros(1, 1);

%% RETRIEVE DATA
% Constants
xpdotp = 1440 / (2*pi);

for i = 1:N
    % Read and divide into lines data
    rawData = splitlines(webread(url{i}));
    S = size(rawData, 1);
    for j = 1:S
        if not(isempty(strtrim(rawData{j})))
            
            if mod(j-1, 3) == 0 % Retrieve the name of the satellite
                ctr = ctr + 1;
                satData.Name{ctr} = strtrim(rawData{j});
                
            elseif mod(j-2, 3) == 0 % Retrieve the first TLE line of the satellite
                line1 = strtrim(rawData{j});
                for jj = 11:16
                    if (line1(jj) == ' ')
                        line1(jj) = '_';
                    end
                end
                
                if (line1(45) ~= ' ')
                    line1(44) = line1(45);
                end
                line1(45) = '.';
                
                if (line1(8) == ' ')
                    line1(8) = 'U';
                end

                if (line1(10) == ' ')
                    line1(10) = '.';
                end

                for jj = 46:50
                    if (line1(jj) == ' ')
                        line1(jj) = '0';
                    end
                end
                if (line1(52) == ' ')
                    line1(52) = '0';
                end
                if (line1(54) ~= ' ')
                    line1(53) = line1(54);
                end
                line1(54) = '.';
                
                if (line1(63) == ' ')
                    line1(63) = '0';
                end

                if ((length(line1) < 68) || (line1(68) == ' '))
                    line1(68) = '0';
                end
                
                satData.satnum(ctr) = str2double(line1(3:7));
                satData.classification{ctr} = line1(8);
                satData.intldesg{ctr} = strtrim(line1(10:17));
                satData.epochyr(ctr) = str2double(line1(19:20));
                satData.epochdays(ctr) = str2double(line1(21:32));
                satData.ndot(ctr) = str2double(line1(34:43));
                satData.nddot(ctr) = str2double(line1(44:50));
                satData.nexp(ctr) = str2double(line1(51:52));
                satData.bstar(ctr) = str2double(line1(53:59));
                satData.ibexp(ctr) = str2double(line1(60:61));
                satData.elnum(ctr) = str2double(line1(65:68));
                
                satData.nddot(ctr) = satData.nddot(ctr) * 10^satData.nexp(ctr);
                satData.bstar(ctr) = satData.bstar(ctr) * 10^satData.ibexp(ctr);
                
                satData.ndot(ctr) = satData.ndot(ctr) / (xpdotp * 1440);         % [rad/min^2]
                satData.nddot(ctr) = satData.nddot(ctr) / (xpdotp * 1440 * 1440);  % [rad/min^3]
                
                year(ctr) = satData.epochyr(ctr) + 2000;
                [month(ctr), day(ctr), hr(ctr), minute(ctr), sec(ctr)] = days2mdh(year(ctr), satData.epochdays(ctr));
                satData.mjd2000satepoch(ctr) = date2mjd2000([year(ctr) month(ctr) day(ctr) hr(ctr) minute(ctr) sec(ctr)]);
                
            elseif mod(j-3, 3) == 0 % Retrieve the second TLE line of the satellite
                line2 = strtrim(rawData{j});
                
                line2(26) = '.';
     
                for jj = 27:33
                    if (line2(jj) == ' ')
                        line2(jj) = '0';
                    end
                end
    
                satData.inclo(ctr) = str2double(line2(8:16));
                satData.nodeo(ctr) = str2double(line2(17:25));
                satData.ecco(ctr) = str2double(line2(26:33));
                satData.argpo(ctr) = str2double(line2(34:42));
                satData.mo(ctr) = str2double(line2(43:51));
                satData.no_kozai(ctr) = str2double(line2(52:63));
                satData.revnum(ctr) = str2double(line2(64:68));
                
                satData.no_kozai(ctr) = satData.no_kozai(ctr) / xpdotp;
                
                satData.inclo(ctr) = satData.inclo(ctr) * pi/180;
                satData.nodeo(ctr) = satData.nodeo(ctr) * pi/180;
                satData.argpo(ctr) = satData.argpo(ctr) * pi/180;
                satData.mo(ctr) = satData.mo(ctr) * pi/180;
            
                
            end
            
            
        end
    end
    
end


%% INITIALIZE VARIABLES FOR SGP4/SDP4 MODEL
[satData] = sgp4Init(satData, satData.mjd2000satepoch);
parsedData = struct2table(satData);
clear('satData');
satData.data = parsedData;
satData.lastUpdate = datestr(datetime('now'));

%% SAVE TLEs DATA
save(strcat(pwd, '\functions\initialize\TLEs.mat'), 'satData')

data

