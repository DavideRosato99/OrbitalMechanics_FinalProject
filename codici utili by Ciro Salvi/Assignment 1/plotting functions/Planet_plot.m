function [x,y,z] = Planet_plot( planet_ID, coordinates, radius_factor )
%Planet_plot function that plot a planet
%
% PROTOTYPE:
%   [x,y,z] = Planet_plot( planet_ID, coordinates, radius_factor )
%
% INPUT:
%   planet_ID[1]      indices of the planets to be plotted
%                        1:   Mercury
%                        2:   Venus
%                        3:   Earth
%                        4:   Mars
%                        5:   Jupiter
%                        6:   Saturn
%                        7:   Uranus
%                        8:   Neptune
%                        9:   Pluto
%   coordinates[1]    position of the planet's centre [km]
%   radius_factor[1]  scaling factor to see the planet
%
% OUTPUT:
%   x, y, z [1]   coordinates of the sphere's points [km]
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1;     % globe transparency level, 1 = opaque, through 0 = invisible

switch(planet_ID)
    case 1
        image_file  = 'Mercury_texture.jpg';
        mean_radius = astroConstants(21); 
    case 2
%         image_file  = 'Venus_texture.jpg';
%         mean_radius = astroConstants(22); 
        image_file = 'VenusAtmo_texture.jpg';
        mean_radius = astroConstants(22) + 1.6*30; % https://solarsystem.nasa.gov/planets/venus/in-depth/
    case 3
        image_file  = 'Earth_texture.jpg';
        mean_radius = astroConstants(23); 
    case 4
        image_file  = 'Mars_texture.jpg';
        mean_radius = astroConstants(24); 
    case 5
        image_file  = 'Jupiter_texture.jpg';
        mean_radius = astroConstants(25); 
    case 6
        image_file  = 'Saturn_texture.jpg';
        mean_radius = astroConstants(26); 
    case 7
        image_file  = 'Uranus_texture.jpg';
        mean_radius = astroConstants(27); 
    case 8
        image_file  = 'Neptune_texture.jpg';
        mean_radius = astroConstants(28); 
    case 9
        image_file  = 'Moon_texture.jpg';
        mean_radius = astroConstants(30); 
    case 10
        image_file  = 'Sun_texture.jpg';
        mean_radius = astroConstants(3); 
end
% Mean spherical planet
if(nargin <= 2)
    radius_factor = 1;end
erad    = radius_factor*mean_radius; % equatorial radius (km)
prad    = radius_factor*mean_radius; % polar radius (km)
if(nargin == 1 )
    coordinates = [ 0 0 0 ];end
%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(coordinates(1), coordinates(2), coordinates(3), erad, erad, prad, npanels);
hold on
globe = surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

%% Texturemap the globe
% Load planet image for texture map
cdata = imread(image_file);
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
end