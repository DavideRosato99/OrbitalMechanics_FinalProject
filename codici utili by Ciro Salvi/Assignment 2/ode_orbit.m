function dy = ode_orbit(~, y, mu)
% 
% PROTOTYPE:
%   dy = ode_orbit_perturbated(~, y, mu,R_e,J_2)
% 
% INPUT
%   y[6x1] state (position(km), velocity(km/s))
%   mu[1] gravitational constant, km^3/s^2
%   
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
    dy = [ y(4); y(5); y(6); (-mu/r^3)*y(1); (-mu/r^3)*y(2); (-mu/r^3)*y(3)];    
end