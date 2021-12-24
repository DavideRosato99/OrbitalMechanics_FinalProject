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
green0 = wrapTo2Pi(omE * (date0mjd2000*24*60*60))   % [rad] Greenwich 0

%% UNPERTURBED
[raUnp, decUnp, lonUnp, latUnp] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'unpert');

%% REPEATING UNPERTURBED
[raUnpRep, decUnpRep, lonUnpRep, latUnpRep] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'unpert', k, m);

%% PERTURBED
[raPer, decPer, lonPer, latPer] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'pert', date0);

%% REPEATING PERTURBED
% [raPerRep, decPerRep, lonPerRep, latPerRep] = GroundTrackCalc(orbIn, periods, green0, muE, omE, 'pert', date0, k, m);

%% SAVE OUTPUT VARIABLES
%%% UNPERTURBED
data.groundtracks.raUnp = raUnp;
data.groundtracks.decUnp = decUnp;
data.groundtracks.lonUnp = lonUnp;
data.groundtracks.latUnp = latUnp;
%%% UNPERTURBED
data.groundtracks.raUnpRep = raUnpRep;
data.groundtracks.decUnpRep = decUnpRep;
data.groundtracks.lonUnpRep = lonUnpRep;
data.groundtracks.latUnpRep = latUnpRep;
%%% UNPERTURBED
data.groundtracks.raPer = raPer;
data.groundtracks.decPer = decPer;
data.groundtracks.lonPer = lonPer;
data.groundtracks.latPer = latPer;
%%% UNPERTURBED
data.groundtracks.raUnpRep = raUnpRep;
data.groundtracks.decUnpRep = decUnpRep;
data.groundtracks.lonUnpRep = lonUnpRep;
data.groundtracks.latUnpRep = latUnpRep;

%% PLOT
if settings.groundtracks.plot
    for i = 1:length(periods)
        %%% UNPERTURBED
        if i == 1
            str = '1 Period';
        else
            days = periods(i)/(24*60*60);
            str = strcat(num2str(days), " days");
        end
        figure('Name', strcat("Unperturbed - ", str), 'NumberTitle', 'off');
        hold on;
        axis equal;
        set(gca,'XTick',[-180:30:180],'XTickMode','manual');
        set(gca,'YTick',[-90:30:90],'YTickMode','manual');
        xlim([-180,180]); ylim([-90,90]);
        
        image_file = 'earth.png';
        cdata      = flip(imread(image_file));
        imagesc([-180,180],[-90, 90],cdata);
        

        
        plot(lonUnp{i}/pi*180,latUnp{i}/pi*180,'.g')
        hold on 
        plot(lonUnp{i}/pi*180,latUnp{i}/pi*180,'.r','Marker','square','MarkerIndices',1);
        plot(lonUnp{i}(1)/pi*180,latUnp{i}(1)/pi*180,'ro','LineWidth',2);
        plot(lonUnp{i}(end)/pi*180,latUnp{i}(end)/pi*180,'rs','LineWidth',2);
        
        
        

    end

end


end