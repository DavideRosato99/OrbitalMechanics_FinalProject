function [errorFB, r_p, h_ga, delta, delta_V_powFB, e_minus, e_plus,...
    a_minus, a_plus] = poweredFlyby(VV_minus, VV_plus, VV_p, h_atm, ID)
% poweredFlyby - Function to obtain the characteristics of a powered flyby
%
% PROTOTYPE
%   [errorFB,r_p,h_ga,delta,delta_V_powFB,e_minus,e_plus,a_minus,a_plus]...
%   =poweredFlyby(VV_minus,VV_plus,VV_p,h_atm,ID)
%
% INPUT:
%   VV_minus         double  [3x1]   s/c velocity on the incoming 
%                                    hyperbola in heliocentric frame [km/s]   
%   VV_plus          double  [3x1]   s/c velocity on the outgoing
%                                    hyperbola in heliocentric frame [km/s] 
%   VV_plus          double  [3x1]   planet velocity in heliocentric 
%                                    frame                           [km/s] 
%   h_atm            double  [1x1]   min altitude of pow. fly-by     [km/s]
%   ID               double  [1x1]   fly-by planet identifier        [-]
%
% OUTPUT:
%   errorFB          double  [1x1]   error in computing the radius of
%                                    pericenter (1: error)           [-]
%   r_p              double  [1x1]   perigree radius                 [km]   
%   h_ga             double  [1x1]   altitude of perigree            [km]
%   Delta            double  [1x1]   turning angle                   [deg]
%   delta_V_powFB    double  [1x1]   cost of the Fly-By manouvre     [km/s]
%   e_minus          double  [1x1]   incoming eccentricity           [-]
%   e_plus           double  [1x1]   outgoing eccentricity           [-]
%   a_minus          double  [1x1]   incoming semi-major axis        [km]
%   a_plus           double  [1x1]   outgoing semi-major axis        [km]
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

warning off

%%
R_planet = astroConstants(20 + ID);
mu = astroConstants(10 + ID);

%% v_minus and v_plus
vv_inf_minus = VV_minus - VV_p;
v_inf_minus = norm(vv_inf_minus);

vv_inf_plus = VV_plus - VV_p;
v_inf_plus = norm(vv_inf_plus);

%% Delta
delta = acos(dot(vv_inf_minus, vv_inf_plus)/(v_inf_minus*v_inf_plus));

%% rp
delta_minus  = @(r_p) 2*asin(1./(1+ r_p * v_inf_minus^2/mu));
delta_plus   = @(r_p) 2*asin(1./(1+ r_p * v_inf_plus^2/mu));
r_p_SolveFun = @(r_p) abs((delta_minus(r_p) + delta_plus(r_p))/2 - delta);

r_p_min = R_planet + h_atm;

options = optimoptions('fsolve','TolFun',1e-14,'Display','off');
r_p = fsolve(r_p_SolveFun, r_p_min, options);  % Periapsis radius of the hyperbolas [km]

if r_p < r_p_min
    
    errorFB = 1;

    delta_V_powFB = NaN;

    h_ga = NaN;

    e_minus = NaN;
    e_plus  = NaN;
    a_minus = NaN;
    a_plus  = NaN;
    
else
    
    errorFB = 0;

    v_p_minus = sqrt(v_inf_minus^2 + 2*mu/r_p);
    v_p_plus  = sqrt(v_inf_plus^2 + 2*mu/r_p);

    delta_V_powFB = abs(v_p_plus - v_p_minus);

    h_ga = r_p - R_planet;

    e_minus = 1/sin(delta_minus(r_p)/2);
    e_plus = 1/sin(delta_plus(r_p)/2);
    a_minus = r_p/(1 - e_minus);
    a_plus = r_p/(1 - e_plus);

end