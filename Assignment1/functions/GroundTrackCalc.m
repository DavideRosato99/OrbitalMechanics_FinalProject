function [ra, dec, lon, lat] = GroundTrackCalc(orbIn, tMaxSpan, green0, mu, om, type, varargin)
% GroundTrack - the function computes the ground tracks of the spacecraft.
%
% PROTOTYPE
%   [ra,dec,lon,lat]=GroundTrack(orbIn,tMaxSpan,green0,mu,om,'unpert') will
%                    perform an unperturbed groundtrack calculation.
% 
%   [ra,dec,lon,lat]=GroundTrack(orbIn,tMaxSpan,green0,mu,om,'pert') will
%                    perform a perturbed groundtrack calculation.
%
%   [ra,dec,lon,lat]=GroundTrack(orbIn,tMaxSpan,green0,mu,om,'unpert', k,
%                    m, date0) will perform an unperturbed repeating
%                    groundtrack calculation.
%
%   [ra,dec,lon,lat]=GroundTrack(orbIn,tMaxSpan,green0,mu,om,'pert', k, m,
%                    date0) will perform a perturbed repeating groundtrack
%                    calculation.
%
% INPUT:
%   orbIn      double  [1x6]   Initial orbit parameter           [s]
%   tMaxSpan   double  [Nx6]   Maximum time for calculation      [s]
%   green0     double  [1x1]   Starting Greenwich                [rad]
%   mu         double  [1x1]   Planet Gravitational parameter    [km^3/s^2]
%   om         double  [1x1]   Planet angular velocity           [rad/s]
%   type       char    [1x1]   Simulation type                   [-]
%
% PROTOTYPE OUTPUT: 
%   ra         double  [1x1]   Rigth ascension                   [rad]
%   dec        double  [1x1]   Declination                       [rad]
%   lon        double  [1x1]   Longitude                         [rad]
%   lat        double  [1x1]   Latitude                          [rad]
%   
% CALLED FUNCTIONS: 

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

%% PRE-CALCULATIONS
N   = length(tMaxSpan);
ra  = cell(N, 1);
dec = cell(N, 1);
lon = cell(N, 1);
lat = cell(N, 1);
deltaT = 10;
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

[rr0, vv0] = par2car(orbIn, mu);
Y0 = [rr0; vv0];

%% CALCULATIONS
for i = 1:N
    t_vec = 0 : deltaT : tMaxSpan(i);
    
    if strcmp(type, 'unpert') && not(isempty(varargin)) % Unperturbed repeating
        k = varargin{1};
        m = varargin{2};
        Tp = 2*pi/(om) * (m/k);
        a_rep = (((Tp/(2*pi))^2)*mu)^(1/3);
        orbIn(1) = a_rep;
        t_vec = 0 : deltaT : Tp;
        
        [rr0, vv0] = par2car(orbIn, mu);
        Y0 = [rr0; vv0];
        
    elseif strcmp(type, 'pert') && not(isempty(varargin))
        
        if size(varargin, 2) == 1 % Perturbed
            date0 = varargin{1};
        else                      % Perturbed repeating
            date0 = varargin{1};
        end
        
    end
    
    %%% Ode
    if strcmp(type, 'unpert')
        [T, Y] = ode113(@ode_2bp, t_vec, Y0, options, mu, 'cart');
    else
        [T, Y] = ode113(@ode_2bp, t_vec, Y0, options, mu, 'cart', datetime(date0));
    end
    
    theta = wrapTo2Pi(green0 + om*T);

    Nt = length(T);

    ra{i}  = zeros(Nt, 1);
    dec{i} = zeros(Nt, 1);
    lat{i} = zeros(Nt, 1);
    lon{i} = zeros(Nt, 1);
    
    for j = 1:Nt
        R = [cos(theta(j)), sin(theta(j)), 0;
            -sin(theta(j)), cos(theta(j)), 0;
                   0              0        1];

        r_local = R*Y(j,1:3)';

        r = norm(Y(j,1:3));

        dec{i}(j) = asin(Y(j,3)/r);

        if Y(j,2) > 0
            ra{i}(j) = acos(Y(j,1)/r/cos(dec{i}(j)));
        else
            ra{i}(j) = 2*pi-acos(Y(j,1)/r/cos(dec{i}(j)));
        end

        lat{i}(j) = asin(r_local(3)/r);

        if r_local(2) > 0
            lon{i}(j) = acos(r_local(1)/r/cos(lat{i}(j)));
        else
            lon{i}(j) = 2*pi-acos(r_local(1)/r/cos(lat{i}(j)));
        end

        if lon{i}(j) > pi
            lon{i}(j) = lon{i}(j) - 2*pi;
        end

    end

    
end



end