function dy = ode_orbit(~, y, mu)
%ode_orbit ODE system
%
% PROTOTYPE:
%   dy = ode_orbit(~, y, mu)
%
% INPUT:
%   t[1]        Time (can be omitted, as the system is autonomous) [T]
%   y[6x1]      Cartesian state of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T]
%   muP [1]     Gravitational parameter of the primary [L^3/T^
%
% OUTPUT:
%   dy [6x1]    Derivative of the state [ L/T^2, L/T^3]
%
% CONTRIBUTORS:
% Student 1
% Student 2
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%   Suhailah Alkhawashke
%   Ciro Salvi
%
% VERSIONS
%   2021-02-11
%

    pos = [y(1) y(2) y(3)];
    r = norm(pos);
    dy = [ y(4); y(5); y(6); (-mu/r^3)*y(1); (-mu/r^3)*y(2); (-mu/r^3)*y(3)];    
end