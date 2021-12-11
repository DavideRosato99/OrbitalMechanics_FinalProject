function [DeltaV, h_GA, r_p, v_p_minus, v_p_plus] = poweredGA(fb_ID, VV_p, VV_minus, VV_plus)
% PROTOTYPE:
%  [DeltaV, h_GA, r_p, v_p_minus, v_p_plus] = poweredGA(fb_ID, VV_p, VV_minus, VV_plus)
% 
% DESCRIPTION:
%   Returns relevant parameters of the Powered Gravity Assist
% 
% INPUT:
%     fb_ID: Integer number identifying the celestial body (< 10)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%     VV_p: planet velocity vector (heliocentric frame)
%     VV_minus:  s/c velocity vector before the flyby(heliocentric frame)
%     VV_plus: s/c velocity vector before the flyby(heliocentric frame)
% 
% OUTPUT:
%     DeltaV: Delta velocity of the Power Gravity Assist
%     h_GA: Altitude of the closest approach with the planet
%     r_p:  Periapsis radius of the hyperbolas 
%     v_p_minus, v_p_plus: velocity before and after the flyby wrt
%     planet velocity

% CALLED FUNCTIONS:
%     astroConstants.m

mu_s = astroConstants(4); % [km^3/s^2]
mu_p = astroConstants(10 + fb_ID);
R_planet = astroConstants(20 + fb_ID); % max radius (i.e. equatorial)
h_atm = 2000;
r_p_min = R_planet + h_atm;


vv_inf_minus = VV_minus - VV_p;
v_inf_minus = norm(vv_inf_minus);

vv_inf_plus = VV_plus - VV_p;
v_inf_plus = norm(vv_inf_plus);

delta = acos(dot(vv_inf_minus, vv_inf_plus)/(v_inf_minus*v_inf_plus));



delta_minus  = @(r_p) 2*asin(1./(1+ r_p * v_inf_minus^2/mu_p));
delta_plus   = @(r_p) 2*asin(1./(1+ r_p * v_inf_plus^2/mu_p));
r_p_SolveFun = @(r_p) (delta_minus(r_p) + delta_plus(r_p))/2 - delta;
options = optimset('TolFun',1e-14,'Display','off');


r_p = fzero(r_p_SolveFun, r_p_min, options);  % Periapsis radius of the hyperbolas [km]

v_p_minus = sqrt(v_inf_minus^2 + 2*mu_p/r_p);
v_p_plus = sqrt(v_inf_plus^2 + 2*mu_p/r_p);

DeltaV = abs(v_p_plus - v_p_minus);


h_GA = r_p - R_planet;


if isnan(r_p) || r_p < r_p_min
    DeltaV = nan;
end


end

