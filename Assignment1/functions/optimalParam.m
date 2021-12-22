function [orbOpt] = optimalParam(data, satData, settings)
% optimalParam - Function to find the optimal RAAN and om to ensure maximum
% distances from other objects orbiting around Earth. The latter are
% computed through TLEs.
%
% PROTOTYPE
%   [Re] = optimalParam(data, satData, settings)
%
% INPUT:
%   data       struct  [1x1]   general data struct                  [-]
%   satData    struct  [1x1]   satellites data from TLEs            [-]
%   settings   struct  [1x1]   settings struct                      [-]
%
% OUTPUT:
%   dx       double [6x1]   state vector derivative                 [-]
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

%% PRE-CALCULATION VARIABLES SET UP
%%% Starting
date0  = data.starting.date;              % [-] Satellite departure date
Torbit = data.starting.Torbit;            % [s] Given orbit period
orbIn  = data.starting.orbIn;             % [-] Given orbit initial parameters
%%% Optimal
nPeriod = data.optimal.nPeriod;           % [-] Number of periods the orbit is propagated for
nPoints = data.optimal.nPoints;           % [-] Number of points for each period at which coordinates are computed
nOM     = data.optimal.nOM;               % [-] Number of RAAN that are calculated to find the optimal value
nom     = data.optimal.nom;               % [-] Number of omega that are calculated to find the optimal value
%%% Constants
muE = data.constants.muE;                 % [km^3/s^2] Planetary constant of the Earth
%%% TLEs
satData = satData.data;                   % [-] TLEs


%% COMPUTING POSITION VECTORS OF TLEs -------------------------------------
% Get all earth orbiting objects positions and velocities (TEME reference frame)
[rTEME, vTEME] = SGP4(satData, date0);
[rECI, ~] = teme2eci(rTEME, 0, date2mjd2000(date0));

%% PARSING TLEs THROUGH THEIR APOGEES AND PERIGEES
N   = size(satData, 1);
apoTLE = zeros(N, 1);
perTLE = zeros(N, 1);

apoSat = orbIn(1) * (1 + orbIn(2));
perSat = orbIn(1) * (1 - orbIn(2));
for i = 1:N
    [orb] = car2par(rTEME(i, :), vTEME(i, :), muE);
    apoTLE(i) = orb(1) * (1 + orb(2));
    perTLE(i) = orb(1) * (1 - orb(2));
end

indexesTLE = apoTLE >= perSat & perTLE <= apoSat;
Ntle = sum(indexesTLE);
indContr = zeros(length(indexesTLE), 1);
ctr = 1;
for i = 1:length(indexesTLE)
    if indexesTLE(i)
        indContr(i) = ctr;
        ctr = ctr + 1;
    end
end


%% PROPAGATION OF SELECTED TLEs OBJECTS ----------------------------------------
[Y, MO, D] = ymd(datetime(date0));
[H, M, S] = hms(datetime(date0));

% Define timespan
timespan = linspace(0, nPeriod*Torbit, nPeriod*nPoints);
RtleTEME = zeros(Ntle, 3, length(timespan));
RtleECI  = zeros(Ntle, 3, length(timespan));

if timespan(end)/(24*60*60) > 30
    warning('Total number of days for analysis exceed 30. Please decrease nPeriod')
end

fprintf('** FIND OPTIMAL RAAN AND OM **\n')
fprintf('Total days for analysis:                  %2.1f\n', timespan(end)/(24*60*60))
fprintf("Propagation of all TLEs objects:            ");
if settings.optimal.parallel
    parforProgress(length(timespan));
    parfor i = 1:length(timespan)
        t = timespan(i);
        date = [Y MO D H M S+t];
        [Yl, MOl, Dl] = ymd(datetime(date));
        [Hl, Ml, Sl] = hms(datetime(date));
        % Get TLEs coordinates using SGP4 propagation algorithm
        [RtleTEME(:, :, i), ~] = SGP4(satData(indexesTLE, :), [Yl MOl Dl Hl Ml Sl]);
        [RtleECI(:, :, i), ~] = teme2eci(RtleTEME(:, :, i), 0, date2mjd2000([Yl MOl Dl Hl Ml Sl]));
        fprintf(repmat('\b', 1, 5));
        parforProgress       
    end
    parforProgress(0);
else
    fprintf(repmat('\b', 1, 2));
    str = fprintf('0%%');
    for i = 1:length(timespan)
        t = timespan(i);
        date = [Y MO D H M S+t];
        [Yl, MOl, Dl] = ymd(datetime(date));
        [Hl, Ml, Sl] = hms(datetime(date));
        % Get TLEs coordinates using SGP4 propagation algorithm
        [RtleTEME(:, :, i), ~] = SGP4(satData(indexesTLE, :), [Yl MOl Dl Hl Ml Sl]);
        [RtleECI(:, :, i), ~] = teme2eci(RtleTEME(:, :, i), 0, date2mjd2000([Yl MOl Dl Hl Ml Sl]));
        fprintf(repmat('\b', 1, str));
        str = fprintf(strcat(num2str(i*100/length(timespan)), '%%'));
    end
end



%% SATELLITE UNPERTURBED ORBIT PROPAGATION --------------------------------
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
OmLoop  = deg2rad(linspace(1, 360, nOM));
omLoop  = deg2rad(linspace(1, 360, nom));
I = length(OmLoop);
J = length(omLoop);
Yorb = zeros(length(timespan), 3, I*J);

% Propagat the initial orbit
orb = [orbIn(1:2), 0, 0, 0, orbIn(6)];
[rr, vv] = par2car(orb, muE);
x0 = [rr; vv];
[T0, Y0] = ode113(@ode_2bp, timespan, x0, options, muE, 'cart');

R_i = [1 0 0; 0 cos(orbIn(3)) sin(orbIn(3)); 0 -sin(orbIn(3)) cos(orbIn(3))];
ctr = 1;

fprintf('Satellite unperturbed orbit propagation:  ');
str = fprintf("0 %%");
for i = 1:I
    Om = OmLoop(i);
    R_OM = [cos(Om) sin(Om) 0; -sin(Om) cos(Om) 0; 0 0 1];
    for j = 1:J
        om = omLoop(j);
        R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
        R = R_om * R_i * R_OM;
        Yorb(:, :, ctr) = (R'*Y0(:, 1:3)')';
        ctr = ctr+1;
    end
    fprintf(repmat('\b', 1, str));
    str = fprintf(strcat(num2str((i)*100/I), '%%'));
end
fprintf('\n');


%% FIND MINIMUM DISTANCES -------------------------------------------------
MINdis        = zeros(I, J);
MINdisDEC     = zeros(I, J);
indWorstTLE   = zeros(I, J);

% Take into account the instant of the closest approach
f = @(x) -(1/(nPeriod*Torbit)^2)*x.^2 + (2/(nPeriod*Torbit))*x;

fprintf('Find minimum distances:                     ');
if settings.optimal.parallel
    parforProgress(I);
    parfor i = 1:I
        for j = 1:J
            minDis = zeros(Ntle, 1);
            ind = zeros(Ntle, 1);
            for m = 1:Ntle
                [minDis(m), ind(m)] = min(vecnorm(squeeze(...
                    Yorb(:, 1:3, (i-1)*J + j))' - ...
                    squeeze(RtleECI(m, :, :))));
            end
            [MINdis(i, j), ind1] = min(minDis);
            MINdisDEC(i, j) = min(MINdis(i, j) * f(T0(ind(ind1))));
            indWorstTLE(i, j) = ind1;
        end
        fprintf(repmat('\b', 1, 5));
        parforProgress
    end  
    parforProgress(0);
else
    str = fprintf("0 %%");
    for i = 1:I
        for j = 1:J
            minDis = zeros(Ntle, 1);
            ind = zeros(Ntle, 1);
            for m = 1:Ntle
                [minDis(m), ind(m)] = min(vecnorm(squeeze(...
                    Yorb(:, 1:3, (i-1)*J + j))' - ...
                    squeeze(RtleECI(m, :, :))));
            end
            [MINdis(i, j), ind1] = min(minDis);
            MINdisDEC(i, j) = min(MINdis(i, j) * f(T0(ind(ind1))));
            indWorstTLE(i, j) = ind1;
        end
        fprintf(repmat('\b', 1, str));
        str = fprintf(strcat(num2str((i)*100/I), '%%'));
    end
end
fprintf('\n');


%% FIND OPTIMAL ORBIT AND PROPAGATE IT
BESTdec = max(max(MINdisDEC));
[i1, i2] = find(MINdisDEC == BESTdec);
OmBest = OmLoop(i1)*180/pi;
omBest = omLoop(i2)*180/pi;
fprintf('OPTIMAL OM and om with parabolic decrement function found!\nOM =   %.2f\nom =   %.2f\n',...
    OmBest, omBest);

% Save optimal orbit
orbOpt = [orbIn(1:3), OmLoop(i1), omLoop(i2), orbIn(6)];

% Propagate orbit
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
timespan = linspace(0, nPeriod*Torbit, nPeriod*1000);
[rr, vv] = par2car(orbOpt, muE);
x0 = [rr; vv];

[~, YorbOpt] = ode113(@ode_2bp, timespan, x0, options, muE, 'cart');


%% FIND WORST TLE AND PROPAGATE IT
worstTLE = satData((indContr == indWorstTLE(i1, i2)), :);
fprintf('WORST OBJECT:     %s\n', string(worstTLE.Name));

% Propagate TLE
[Y, MO, D] = ymd(datetime(date0));
[H, M, S]  = hms(datetime(date0));
RworstTLE_TEME = zeros(length(timespan), 3);
RworstTLE_ECI = zeros(length(timespan), 3);

for i = 1:length(timespan)
    t = timespan(i);
    date = [Y MO D H M S+t];
    [Yl, MOl, Dl] = ymd(datetime(date));
    [Hl, Ml, Sl] = hms(datetime(date));
    % Get TLEs coordinates using SGP4 propagation algorithm
    [RworstTLE_TEME(i, :), ~] = SGP4(satData((indContr == indWorstTLE(i1, i2)), :),...
        [Yl MOl Dl Hl Ml Sl]);
    [RworstTLE_ECI(i, :), ~] = teme2eci(RworstTLE_TEME(i, :), 0, ...
        date2mjd2000([Yl MOl Dl Hl Ml Sl]));
end



%% PLOT
if settings.optimal.plot
    %%% LOAD EARTH TEXTURE
    Re = astroConstants(23);
    [Xe, Ye, Ze] = sphere(100);
    Xe = Re*Xe; Ye = Re*Ye; Ze = Re*Ze;
    earthPNG = imread('earth.png');
 
    %%% TOTAL TLEs AROUND EARTH
    figure('Name', 'Objects orbiting the earth', 'NumberTitle', 'off');
    surf(Xe, Ye, Ze, 'CData', flipud(earthPNG), 'FaceColor', 'texture', 'edgecolor', 'none');
    hold on; grid on; axis vis3d
    plot3(rECI(:, 1), rECI(:, 2), rECI(:, 3), 'o', 'MarkerSize', 0.7)
    xlim([-5e4 5e4]); ylim([-5e4 5e4]); zlim([-5e4 5e4])
    view(45, 20)
    
    %%% TOTAL ORBITS COMPUTED
    figure('Name', 'Total orbits computed', 'NumberTitle', 'off');
    surf(Xe, Ye, Ze, 'CData', flipud(earthPNG), 'FaceColor', 'texture', 'edgecolor', 'none');
    hold on; grid on; axis equal;
    color = parula(I*J);
    for i = 1 : I*J
        plot3(Yorb(:, 1, i), Yorb(:, 2, i), Yorb(:, 3, i), 'color', color(i, :), 'LineWidth', 0.1);
    end
    
    %%% PARABOLIC DECREMENT
    figure('Name', 'Parabolic decrement', 'NumberTitle', 'off');
    x = linspace(0, nPeriod*Torbit, 2000);
    plot(x/(24*60*60), f(x)); grid on;
    xlabel('days [-]'); ylabel('f(x) [-]');
    title('Parabolic decrement', 'interpreter', 'latex');
    
    %%% DISTANCES
    figure('Name','Satellite distance from TLEs','NumberTitle','off');
    [OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
    surf(OMM, OM, MINdis'); hold on;
    shading interp
    xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$min(|\textbf{r} - \textbf{r$_{TLE}$}|)$')
    BEST = max(max(MINdis));
    [i1, i2] = find(MINdis == BEST);
    plot3(OmLoop(i1(1))*180/pi, omLoop(i2(1))*180/pi, BEST, 'ro')
    legend('Distances', 'Optimal', 'interpreter', 'latex')
    
    %%% DISTANCES WITH DECREMENT FUNCTION
    figure('Name','Satellite distance from TLEs with decrement function','NumberTitle','off');
    [OMM, OM] = meshgrid(rad2deg(OmLoop), rad2deg(omLoop));
    surf(OMM, OM, MINdisDEC'); hold on;
    shading interp
    xlabel('$\Omega [deg]$'); ylabel('$\omega [deg]$'); zlabel('$f(\Omega,\omega)$')
    BESTdec = max(max(MINdisDEC));
    [i1, i2] = find(MINdisDEC == BESTdec);
    plot3(OmLoop(i1(1))*180/pi, omLoop(i2(1))*180/pi, BESTdec, 'ro')
    legend('Distances', 'Optimal', 'interpreter', 'latex')
    
    %%% WORST TLE
    figure('Name','Satellite distance from TLEs with decrement function','NumberTitle','off');
    surf(Xe, Ye, Ze, 'CData', flipud(earthPNG), 'FaceColor', 'texture', 'edgecolor', 'none');
    hold on; grid on; axis equal;
    plot3(YorbOpt(:, 1), YorbOpt(:, 2), YorbOpt(:, 3), 'LineWidth', 0.5);
    plot3(RworstTLE_ECI(:, 1), RworstTLE_ECI(:, 2), RworstTLE_ECI(:, 3), 'LineWidth', 0.5);
    sat = plot3(YorbOpt(1, 1), YorbOpt(1, 2), YorbOpt(1, 3), 'bo', 'MarkerSize', 4);
    TLE = plot3(RworstTLE_ECI(1, 1), RworstTLE_ECI(1, 2), RworstTLE_ECI(1, 3), 'ro', 'MarkerSize', 4);
    dis = plot3([YorbOpt(1, 1) RworstTLE_ECI(1, 1)], [YorbOpt(1, 2) RworstTLE_ECI(1, 2)], ...
        [YorbOpt(1, 3) RworstTLE_ECI(1, 3)], 'k-');
    for i = 1:size(YorbOpt, 1)
        delete(sat)
        delete(TLE)
        delete(dis)
        sat = plot3(YorbOpt(i, 1), YorbOpt(i, 2), YorbOpt(i, 3), 'bo', 'MarkerSize', 4);
        TLE = plot3(RworstTLE_ECI(i, 1), RworstTLE_ECI(i, 2), RworstTLE_ECI(i, 3), 'ro', 'MarkerSize', 4);
        dis = plot3([YorbOpt(i, 1) RworstTLE_ECI(i, 1)], [YorbOpt(i, 2) RworstTLE_ECI(i, 2)], ...
        [YorbOpt(i, 3) RworstTLE_ECI(i, 3)], 'k-');
        norm(YorbOpt(i, 1:3) - RworstTLE_ECI(i, :))
        drawnow limitrate
    end
    
    
end





