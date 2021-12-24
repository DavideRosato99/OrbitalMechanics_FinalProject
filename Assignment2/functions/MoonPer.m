function [apMG_vect] = MoonPer(r_moon_vect,r_vect)
%
% Moonper - the function gives the perturbed acceleration due to lunar     
%           gravity influence in cartesian coordinates.
%
% INPUT:
%	r_moon_vect   double  [3x1]    Moon position vector              [km]
%   r_vect        double  [3x1]    Space-craft position vector       [km]
%
% OUTPUT:
%   [apMG_vect]   double  [3x1]    Perturbed acceleration vector [km/s^2]
%
% CALLED FUNCTIONS: 
%   astroConstants
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

mu_moon = astroConstants(20);  % [Km^3/s^2]
r_ms_vect = r_moon_vect - r_vect;

r_ms = norm(r_ms_vect);

apMG_vect = mu_moon*(r_ms_vect/r_ms^3);

end