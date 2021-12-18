function [] = GroundTrack(Tfin, orb, date0, type, varargin)
% GroundTrack - the function computes the ground track of a spacecraft.
%
% PROTOTYPE
%   [ra,dec,lon,lat]=GroundTrack(T,Y,green0,mu)
%
% INPUT:
%   T        double  [1x1]   Time vector                          [s]
%   Y        double  [Nx6]   Position Vector (ECI)               [km]
%   green0   double  [1x1]   Starting Greenwich                 [rad]
%
% PROTOTYPE OUTPUT: 
%   ra       double  [1x1]   rigth ascension                       []
%   dec      double  [1x1]   declination                           []
%   lon      double  [1x1]   longitude                             []
%   lat      double  [1x1]   latitude                              []
%   
% CALLED FUNCTIONS: 
%   astroConstants 
%   par2car
%   ode_2bp
%
% NOTE: 
%   .Time vector must be the ode solution time vector, starting from 0;
%   .If the function is used without output it gives only the plot of the
%    groundtrack.
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

%%
muE = astroConstants(13);
deltaT = 10;
t_vec = 0:deltaT:Tfin;

om_E = (15.04*pi/180)/(60*60);

if not(isempty(varargin))
    m = varargin{1};
    k = varargin{2};
    Te = 2*pi/(om_E) * (m/k);
    a_rep = (((Te/(2*pi))^2)*muE)^(1/3);
    orb(1) = a_rep;
    t_vec = 0:deltaT:Te;
end

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[rr0, vv0] = par2car(orb, muE);
Y0 = [rr0; vv0];
if strcmp(type, 'unpert')
    [T, Y] = ode113(@ode_2bp, t_vec, Y0, options, muE);
else
    [T, Y] = ode113(@ode_2bp, t_vec, Y0, options, muE, date0);
end

deltaSec = hms(date0 - datetime([2000 1 1 12 0 0]));
green0 = wrapTo2Pi((2*pi/(24*60*60)) * deltaSec);
theta = wrapTo2Pi(green0 + om_E*T);

N = length(T);

ra = zeros(N, 1);
dec = zeros(N, 1);
lat = zeros(N, 1);
lon = zeros(N, 1);

for i = 1:N
    R = [cos(theta(i)), sin(theta(i)), 0;
        -sin(theta(i)), cos(theta(i)), 0;
               0              0        1];
           
    r_local = R*Y(i,1:3)';
    
    r = norm(Y(i,1:3));
    
    dec(i) = asin(Y(i,3)/r);
    
    if Y(i,2) > 0
        ra(i) = acos(Y(i,1)/r/cos(dec(i)));
    else
        ra(i) = 2*pi-acos(Y(i,1)/r/cos(dec(i)));
    end
    
    lat(i) = asin(r_local(3)/r);
    
    if r_local(2) > 0
        lon(i) = acos(r_local(1)/r/cos(lat(i)));
    else
        lon(i) = 2*pi-acos(r_local(1)/r/cos(lat(i)));
    end
    
    if lon(i) > pi
        lon(i) = lon(i) - 2*pi;
    end
    
end

figure
hold on;
axis equal;
set(gca,'XTick',[-180:30:180],'XTickMode','manual');
set(gca,'YTick',[-90:30:90],'YTickMode','manual');
xlim([-180,180]); ylim([-90,90]);
        

image_file = 'earth.png';
cdata      = flip(imread(image_file));
imagesc([-180,180],[-90, 90],cdata);

plot(lon/pi*180,lat/pi*180,'.g'); hold on;
s = plot(lon(1)/pi*180,lat(1)/pi*180, '.b', 'MarkerSize', 20);
e = plot(lon(end)/pi*180,lat(end)/pi*180, '.r', 'MarkerSize', 20);
legend([s, e], {'start', 'end'});




end