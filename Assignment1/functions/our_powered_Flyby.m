function [r_p, h_ga, Delta, delta_V_poweredFB,e_minus,e_plus,a_minus,a_plus] = our_powered_Flyby(VV_minus,VV_plus)
mu_p = astroConstants(13); % [km^3/s^2]
mu_s = astroConstants(4); % [km^3/s^2]
AU = astroConstants(2); % [km]
R_planet = astroConstants(23); % max radius (i.e. equatorial)
h_atm = 200; % we still have to decide this

V_p = sqrt(mu_s/AU); % we assume a circualr Earth orbit
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
