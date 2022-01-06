function [delta_t]= FlyBytof(e_minus, a_minus,e_plus,a_plus, ID)
% FlyBytof - Function to compute the time of flight of the Fly By given the
% characteristics of the incoming and outgoing hyperbolas
%
% PROTOTYPE
%   [delta_t]= FlyBytof(e_minus, a_minus,e_plus,a_plus, ID)
%
% INPUT:
%   e_minus   double [1x1]   eccentricity of the incoming hyperbola    [-]
%   a_minus   double [1x1]   semi major axis of the incoming hyperbola [km]
%   e_plus    double [1x1]   eccentricity of the outgoing hyperbola    [-]
%   a_plus    double [1x1]   semi major axis of the outgoing hyperbola [km]
%   date      double [1x1]   date
%   ID        double [1x1]   planet ID linked to astroConstants
%
% OUTPUT:
%   delta_t   double [1x1]   time of flight of the Fly By              [s]
%
% CALLED FUNCTIONS: astroConstants, uplanet, par2car
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
rsoi = astroConstants(2)*((astroConstants(10+ID)/astroConstants(4)))^(2/5);

p_minus = a_minus*(1-e_minus^2);
p_plus = a_plus*(1-e_plus^2);


theta_minus = 2*pi - acos(1/e_minus*((p_minus/rsoi) -1));
theta_plus = acos(1/e_plus*((p_plus/rsoi) -1));

F_minus = 2* atanh(sqrt((e_minus-1)/(e_minus+1))*tan(theta_minus/2));
F_plus = 2* atanh(sqrt((e_plus-1)/(e_plus+1))*tan(theta_plus/2));

n_minus = sqrt(astroConstants(10+ID)/a_minus^3);
n_plus = sqrt(astroConstants(10+ID)/a_plus^3);

deltat_minus = 1/n_minus * (F_minus-e_minus*sinh(F_minus));
deltat_plus = 1/n_plus * (F_plus-e_plus*sinh(F_plus));

delta_t = abs(deltat_minus) + abs(deltat_plus);
delta_t_days = delta_t/(3600*24);

fprintf('Total days in SOI: %g\n\n', delta_t_days);
