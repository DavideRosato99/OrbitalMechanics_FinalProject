function [x,y,z] =  SOI_plot(radius, alpha, coordinates)
%SOI_plot plot the Sphere of Influence around a planet
%
% PROTOTYPE:
%   [x,y,z] =  SOI_plot(radius, alpha, coordinates)
%
% INPUT:
%   radius[1]           radius of the SOI [km] 
%   alpha[1]            sphere transparency level
%   coordinates[3]      position of the planet's centre [km]
%
% OUTPUT:
%   x, y, z [1]         coordinates of the sphere's points [km]
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
if ( nargin <= 2 )
    coordinates = [ 0 0 0 ];end
if ( nargin == 1 ) %sphere transparency level, 1 = opaque, through 0 = invisible
    alpha   = 0.5;end
%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(coordinates(1), coordinates(2), coordinates(3), radius, radius, radius, npanels);
hold on
globe = surf(x, y, z, 'FaceColor', 'none');
set(globe, 'FaceColor', '#EDB120', 'FaceAlpha', alpha, 'EdgeColor', 'none');
end