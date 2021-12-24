function [data] = realData(data, satData, settings)
% realData - Function to 
%
% PROTOTYPE
%   [data] = realData(data, satData, settings)
%
% INPUT:
%   data       struct  [1x1]   general data struct                  [-]
%   satData    struct  [1x1]   satellites data from TLEs            [-]
%   settings   struct  [1x1]   settings struct                      [-]
%
% OUTPUT:
%   data       double  [1x6]   Optimal orbit keplerian parameters   [-]
%
% CALLED FUNCTIONS: date2mjd2000, SGP4, teme2eci, car2par, par2car,
%                   parfor_progress
%
% NOTE: 
%       - if settings.optimal.parallel is set to TRUE, a parallel computing
%         will be perfomed;
%       - if settings.optimal.plot is set to TRUE, plots will be computed
%         and displayed;
%       - if settings.optimal.movie is set to TRUE, movies will be created
%         and saved in folder ..\functions\movies.
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------
muE     = data.constants.muE;
satData = satData.data;
date0   = data.starting.date;
Tmax    = data.realData.Tmax;
nPoints = data.realData.nPoints;


if not(isnan(data.starting.OM) && isnan(data.starting.om))
    orbIn  = data.starting.orbIn;         % [-] Given orbit initial parameters
end
%%% Optimal
if isnan(data.starting.OM) && isnan(data.starting.om)
    orbIn  = data.optimal.orbIn;          % [-] Given orbit initial parameters
end

%% SELECT TLE
[rTEME, vTEME] = SGP4(satData, date0);
[rECI, vECI] = teme2eci(rTEME, vTEME, date2mjd2000(date0));

N   = size(satData, 1);
val = zeros(N, 1);
valMin = 10e9;

for i = 1:N
    [orb] = car2par(rECI(i, :), vECI(i, :), muE);
    val(i) = abs((orbIn(1)-orb(1))/(orbIn(1))) + 100*abs((orbIn(2)-orb(2))/(orbIn(1)));    
   
    if val(i) < valMin
       valMin = val(i);
       index = i;
    end
    
end


fprintf(strcat('Best satellite:', string(satData.Name{index})));
format long
[rTEME, vTEME] = SGP4(satData(index,:), date0);
[rECI, vECI] = teme2eci(rTEME, vTEME, date2mjd2000(date0));
[orb1] = car2par(rECI, vECI, muE);



%% PROPAGATING



timespan = linspace(0,Tmax,nPoints);
[Y, MO, D] = ymd(datetime(date0));
[H, M, S]  = hms(datetime(date0));
R_TLE_TEME = zeros(length(timespan), 3);
R_TLE_ECI = zeros(length(timespan), 3);
orb = zeros(length(timespan), 6);
for i = 1:length(timespan)
    t = timespan(i);
    date = [Y MO D H M S+t];
    [Yl, MOl, Dl] = ymd(datetime(date));
    [Hl, Ml, Sl] = hms(datetime(date));
    % Get TLEs coordinates using SGP4 propagation algorithm
    [R_TLE_TEME(i, :), V_TLE_TEME(i, :)] = SGP4(satData(index, :),...
        [Yl MOl Dl Hl Ml Sl]);
    [R_TLE_ECI(i, :), V_TLE_ECI(i, :)] = teme2eci(R_TLE_TEME(i, :), V_TLE_TEME(i, :), ...
        date2mjd2000([Yl MOl Dl Hl Ml Sl]));

    orb(i, :) = car2par(R_TLE_ECI(i, :), V_TLE_ECI(i, :), muE);


end


options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[TGauss, YGauss] = ode113(@ode_2bp, timespan, orb1, options, muE, 'gauss', datetime(date0));

%% PLOT

    %%% SEMI-MAJOR AXIS
    figure('Name', 'Semi-major axis evolution', 'NumberTitle', 'off');
    plot(timespan, orb(:, 1)); grid on; hold on;
    plot(timespan, YGauss(:, 1));
    xlabel('Time [years]'); ylabel('a [km]'); title('Semi-major axis evolution')
    legend('Real Data', 'Gauss', 'interpreter', 'latex')
    
    %%% ECCENTRICITY
    figure('Name', 'Eccentricity evolution', 'NumberTitle', 'off');
    plot(timespan, orb(:, 2)); grid on; hold on;
    plot(timespan, YGauss(:, 2));
    xlabel('Time [years]'); ylabel('e [-]'); title('Eccentricity evolution')
    legend('Real Data', 'Gauss', 'interpreter', 'latex')
    
    %%% INCLINATION
    figure('Name', 'Inclination evolution', 'NumberTitle', 'off');
    plot(timespan, rad2deg(orb(:, 3))); grid on; hold on;
    plot(timespan, rad2deg(YGauss(:, 3)));
    xlabel('Time [years]'); ylabel('i [deg]'); title('Inclination evolution')
    legend('Real Data', 'Gauss', 'interpreter', 'latex')
    
    %%% RAAN
    figure('Name', 'RAAN evolution', 'NumberTitle', 'off');
    plot(timespan, rad2deg(orb(:, 4))); grid on; hold on;
    plot(timespan, rad2deg(YGauss(:, 4)));
    xlabel('Time [years]'); ylabel('$\Omega$ [deg]'); title('RAAN evolution')
    legend('Real Data', 'Gauss', 'interpreter', 'latex')
    
    %%% ARGUMENT OF PERICENTER
    figure('Name', 'Argument of pericenter evolution', 'NumberTitle', 'off');
    plot(timespan, rad2deg(orb(:, 5))); grid on; hold on;
    plot(timespan, rad2deg(YGauss(:, 5)));
    xlabel('Time [years]'); ylabel('$\omega$ [deg]'); title('Argument of pericenter evolution')
    legend('Real Data', 'Gauss', 'interpreter', 'latex')


































