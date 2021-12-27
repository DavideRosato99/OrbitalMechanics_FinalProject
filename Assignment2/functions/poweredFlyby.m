function [r_p, h_ga, Delta, delta_V_poweredFB,e_minus,e_plus,a_minus,a_plus] = poweredFlyby(data,settings)
% poweredFlyby - Function to obtain the characteristics of a powered flyby
%
% PROTOTYPE
%   [data] = poweredFlyby(data,settings)
%
% INPUT:
%   data               struct  [1x1]   general data struct           [-]
%   settings           struct  [1x1]   settings struct               [-]
%
% OUTPUT:
%   r_p                double  [1x1]   perigree radius               [km]   
%   h_ga               double  [1x1]   altitude of perigree          [km]
%   Delta              double  [1x1]   turning angle                 [deg]
%   delta_V_poweredFB  double  [1x1]   cost of the Fly-By            [km/s]
%   e_minus            double  [1x1]   incoming eccentricity         [-]
%   e_plus             double  [1x1]   outgoing eccentricity         [-]
%   a_minus            double  [1x1]   incoming semi-major axis      [km]
%   a_plus             double  [1x1]   incoming semi-major axis      [km]
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


mu_p = data.constants.muE;      % [km^3/s^2]
mu_s = data.constants.muS;      % [km^3/s^2]
AU   = data.constants.AU;       % [km]
R_planet = data.constants.Re;   % max radius (i.e. equatorial)
h_atm = data.flyby.hAtm;        % we still have to decide this
VV_minus = data.flyby.VVminus;  % [km^2]
VV_plus = data.flyby.VVplus;    % [km^2]


V_p = sqrt(mu_s/AU); % Assuming circular Earth orbit
VV_p = V_p * [1 0 0]';


vv_inf_minus = VV_minus - VV_p;
v_inf_minus = norm(vv_inf_minus);

vv_inf_plus = VV_plus - VV_p;
v_inf_plus = norm(vv_inf_plus);

delta = acos(dot(vv_inf_minus, vv_inf_plus)/(v_inf_minus*v_inf_plus));
Delta=rad2deg(delta);


delta_minus  = @(r_p) 2*asin(1./(1+ r_p * v_inf_minus^2/mu_p));
delta_plus   = @(r_p) 2*asin(1./(1+ r_p * v_inf_plus^2/mu_p));
r_p_SolveFun = @(r_p) (delta_minus(r_p) + delta_plus(r_p))/2 - delta;
r_p_min = R_planet + h_atm;
options = optimoptions('fsolve','TolFun',1e-14,'Display','off');

r_p = fsolve(r_p_SolveFun, r_p_min, options);  % Periapsis radius of the hyperbolas [km]

v_p_minus = sqrt(v_inf_minus^2 + 2*mu_p/r_p);
v_p_plus = sqrt(v_inf_plus^2 + 2*mu_p/r_p);

delta_V_poweredFB = abs(v_p_plus - v_p_minus);

h_ga = r_p - R_planet;

e_minus = 1/sin(delta_minus(r_p)/2)
e_plus = 1/sin(delta_plus(r_p)/2)
a_minus = r_p/(1-e_minus); % note that for an hyperbola: a<0
a_plus = r_p/(1-e_plus);

end