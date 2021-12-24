function data = GroundTrack(data, settings)
% GroundTrack - the function computes the ground tracks of the spacecraft.
%
% PROTOTYPE
%   [aRep,aRepPert]=GroundTrack(data,settings)
%
% INPUT:
%   data       struct  [1x1]   general data struct                  [-]
%   settings   struct  [1x1]   settings struct                      [-]
%
% PROTOTYPE OUTPUT: 
%   aRep       double  [1x1]   Semi-major axis for unperturbed rep.
%                              groundtrack                          [km]
%   aRepPert   double  [1x1]   Semi-major axis for perturbed rep.
%                              groundtrack                          [km]
%   
% CALLED FUNCTIONS: 
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
if not(isnan(data.starting.OM) && isnan(data.starting.om))
    orbIn  = data.starting.orbIn;         % [-] Given orbit initial parameters
end
%%% Optimal
if isnan(data.starting.OM) && isnan(data.starting.om)
    orbIn  = data.optimal.orbIn;          % [-] Given orbit initial parameters
end
%%% Constants
muE = data.constants.muE;                 % [km^3/s^2] Planetary constant of the Earth
%%% Groundtracks
periods = data.groundtracks.periods;      % [s] Periods for which the groundtracks will be displayed
k = data.groundtracks.k;                  % [-] Number of periods of the Earth
m = data.groundtracks.m;                  % [-] Number of periods of the satellite


%% GREENWICH 0
omE = (15.04*pi/180)/(60*60);                        % [rad/s] Earth rotational velocity
date0mjd2000 = date2mjd2000(date0);
green0 = wrapTo2Pi(omE * (date0mjd2000*24*60*60));   % [rad] Greenwich 0

%% UNPERTURBED
fprintf('** GROUNDTRACKS **\n')
fprintf("Unperturbed groundtrack:                    ");
str = fprintf('0%%\n');
[raUnp, decUnp, lonUnp, latUnp, ~, ~] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'unpert');
fprintf(repmat('\b', 1, str));
fprintf('100%%\n');

%% PERTURBED
fprintf("Perturbed groundtrack:                      ");
str = fprintf('0%%\n');
[raPer, decPer, lonPer, latPer, ~, ~] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'pert', date0);
fprintf(repmat('\b', 1, str));
fprintf('100%%\n');

%% REPEATING UNPERTURBED
fprintf("Repeating unperturbed groundtrack:          ");
str = fprintf('0%%\n');
orbIn
[raUnpRep, decUnpRep, lonUnpRep, latUnpRep, aUnpRep, tUnpRep] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'unpert', k, m);
fprintf(repmat('\b', 1, str));
fprintf('100%%\n');
fprintf("Unperturbed repeating semi-major axis:  %g [km]\n", aUnpRep)
fprintf("Unperturbed repeating period:            %g [s]\n", tUnpRep)


%% REPEATING PERTURBED
fprintf("Repeating perturbed groundtrack:            ");
str = fprintf('0%%\n');
[raPerRep, decPerRep, lonPerRep, latPerRep, aPerRep, tPerRep] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'pert', date0, k, m);
fprintf(repmat('\b', 1, str));
fprintf('100%%\n');
fprintf("Perturbed repeating semi-major axis:    %g [km]\n", aPerRep)
fprintf("Perturbed repeating period:              %g [s]\n", tPerRep)

%% SAVE OUTPUT VARIABLES
%%% UNPERTURBED
data.groundtracks.raUnp  = raUnp;
data.groundtracks.decUnp = decUnp;
data.groundtracks.lonUnp = lonUnp;
data.groundtracks.latUnp = latUnp;
%%% UNPERTURBED
data.groundtracks.raUnpRep  = raUnpRep;
data.groundtracks.decUnpRep = decUnpRep;
data.groundtracks.lonUnpRep = lonUnpRep;
data.groundtracks.latUnpRep = latUnpRep;
%%% UNPERTURBED
data.groundtracks.raPer  = raPer;
data.groundtracks.decPer = decPer;
data.groundtracks.lonPer = lonPer;
data.groundtracks.latPer = latPer;
%%% UNPERTURBED
data.groundtracks.raUnpRep  = raUnpRep;
data.groundtracks.decUnpRep = decUnpRep;
data.groundtracks.lonUnpRep = lonUnpRep;
data.groundtracks.latUnpRep = latUnpRep;

%% PLOT
if settings.groundtracks.plot
    for i = 1:length(periods)
        
        %%% UNPERTURBED
        warning off
        if i == 1
            str = '1 Period';
        else
            days = periods(i)/(24*60*60);
            str = strcat(num2str(days), " days");
        end
        figure('Name', strcat("Unperturbed - ", str), 'NumberTitle', 'off');
        title(strcat("Unperturbed - ", str))
        hold on;
        axis equal;
        set(gca,'XTick', [-180:30:180], 'XTickMode','manual');
        set(gca,'YTick', [-90:30:90], 'YTickMode', 'manual');
        xlim([-180,180]); ylim([-90,90]);
        
        image_file = 'earth.png';
        cdata      = flip(imread(image_file));
        imagesc([-180,180],[-90, 90],cdata);
        
        plot(lonUnp{i}/pi*180, latUnp{i}/pi*180, '.g', 'MarkerSize', 0.5)
        hold on
        s = plot(lonUnp{i}(1)/pi*180, latUnp{i}(1)/pi*180, 'ro', 'LineWidth', 2);
        e = plot(lonUnp{i}(end)/pi*180, latUnp{i}(end)/pi*180, 'rs', 'LineWidth', 2);
        legend([s, e], {'Start', 'End', 'interpreter', 'latex'})
        xlabel('longitude [deg]'); ylabel('latitude [deg]')
        
        
        %%% PERTURBED
        warning off
        if i == 1
            str = '1 Period';
        else
            days = periods(i)/(24*60*60);
            str = strcat(num2str(days), " days");
        end
        figure('Name', strcat("Perturbed - ", str), 'NumberTitle', 'off');
        title(strcat("Perturbed - ", str))
        hold on;
        axis equal;
        set(gca, 'XTick', [-180:30:180], 'XTickMode', 'manual');
        set(gca, 'YTick', [-90:30:90], 'YTickMode', 'manual');
        xlim([-180,180]); ylim([-90,90]);
        
        image_file = 'earth.png';
        cdata      = flip(imread(image_file));
        imagesc([-180,180], [-90, 90], cdata);
        
        plot(lonPer{i}/pi*180, latPer{i}/pi*180, '.g', 'MarkerSize', 0.5)
        hold on
        s = plot(lonPer{i}(1)/pi*180, latPer{i}(1)/pi*180, 'ro', 'LineWidth', 2);
        e = plot(lonPer{i}(end)/pi*180, latPer{i}(end)/pi*180, 'rs', 'LineWidth', 2);
        legend([s, e], {'Start', 'End', 'interpreter', 'latex'})
        xlabel('longitude [deg]'); ylabel('latitude [deg]')
        
    end
    
    %%% UNPERTURBED REPEATING
    warning off
    if i == 1
        str = '1 Period';
    else
        days = periods(i)/(24*60*60);
        str = strcat(num2str(days), " days");
    end
    figure('Name', strcat("Unperturbed Repeating - ", str), 'NumberTitle', 'off');
    title(strcat("Unperturbed Repeating - ", str))
    hold on;
    axis equal;
    set(gca,'XTick', [-180:30:180], 'XTickMode','manual');
    set(gca,'YTick', [-90:30:90], 'YTickMode', 'manual');
    xlim([-180,180]); ylim([-90,90]);

    image_file = 'earth.png';
    cdata      = flip(imread(image_file));
    imagesc([-180,180],[-90, 90],cdata);

    plot(lonUnpRep{i}/pi*180, latUnpRep{i}/pi*180, '.g', 'MarkerSize', 0.5)
    hold on
    s = plot(lonUnpRep{i}(1)/pi*180, latUnpRep{i}(1)/pi*180, 'ro', 'LineWidth', 2);
    e = plot(lonUnpRep{i}(end)/pi*180, latUnpRep{i}(end)/pi*180, 'rs', 'LineWidth', 2);
    legend([s, e], {'Start', 'End', 'interpreter', 'latex'})
    xlabel('longitude [deg]'); ylabel('latitude [deg]')
    
    %%% PERTURBED REPEATING
    warning off
    if i == 1
        str = '1 Period';
    else
        days = periods(i)/(24*60*60);
        str = strcat(num2str(days), " days");
    end
    figure('Name', strcat("Perturbed Repeating - ", str), 'NumberTitle', 'off');
    title(strcat("Perturbed Repeating - ", str))
    hold on;
    axis equal;
    set(gca,'XTick', [-180:30:180], 'XTickMode','manual');
    set(gca,'YTick', [-90:30:90], 'YTickMode', 'manual');
    xlim([-180,180]); ylim([-90,90]);

    image_file = 'earth.png';
    cdata      = flip(imread(image_file));
    imagesc([-180,180],[-90, 90],cdata);

    plot(lonPerRep{i}/pi*180, latPerRep{i}/pi*180, '.g', 'MarkerSize', 0.5)
    hold on
    s = plot(lonPerRep{i}(1)/pi*180, latPerRep{i}(1)/pi*180, 'ro', 'LineWidth', 2);
    e = plot(lonPerRep{i}(end)/pi*180, latPerRep{i}(end)/pi*180, 'rs', 'LineWidth', 2);
    legend([s, e], {'Start', 'End', 'interpreter', 'latex'})
    xlabel('longitude [deg]'); ylabel('latitude [deg]')

end


end