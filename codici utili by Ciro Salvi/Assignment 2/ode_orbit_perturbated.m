function dy = ode_orbit_perturbated(~, y, mu,R_e,J_2)

% 
% PROTOTYPE:
%   dy = ode_orbit_perturbated(~, y, mu,R_e,J_2)
% 
% INPUT
%    y[6x1] state (position(km), velocity(km/s))
%   mu[1] gravitational constant, km^3/s^2
%   R_e radius of earth, km
%   J_2 constant of harmonic zonal perturbation
% OUTPUT
% 	dy derivative of the state

% 
% CONTRIBUTORS
%    Alessandro Staffolani
%    Ciro Salvi
% 
% VERSIONS
% 2020-02-11
    pos = [y(1) y(2) y(3)];
    r = norm(pos);
    a_i = (3/2)*((J_2*mu*R_e^2)/r^4)*(y(1)/r)*(5*((y(3)^2)/r^2)-1);
    a_j = (3/2)*((J_2*mu*R_e^2)/r^4)*(y(2)/r)*(5*((y(3)^2)/r^2)-1);
    a_k = (3/2)*((J_2*mu*R_e^2)/r^4)*(y(3)/r)*(5*((y(3)^2)/r^2)-3);
    dy = [ y(4); y(5); y(6); (-mu/r^3)*y(1)+a_i; (-mu/r^3)*y(2)+a_j; (-mu/r^3)*y(3)+a_k];    
end